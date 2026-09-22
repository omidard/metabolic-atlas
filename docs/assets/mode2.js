// Simulate stage: flux feasibility of found pathways for a chosen GEM on an
// editable medium. Owns the media panel, the substrate-as-carbon-source swap,
// and the per-pathway feasibility / pFBA / FVA runs. All solving happens in the
// GLPK worker; every result carries its denominator; absent stays absent.
// The GEM and medium choices are read from and written to the session context.

import { loadGem, loadMedia, fmt, onEditsChanged } from './data.js';
import { getContext, setContext, onContext } from './context.js';
import { getGLPK, buildLP, boundType, maxGrowth, productTarget, pathwayFeasibility, pathwayPFBA, pathwayFVA, stepConstraints, stepRealizations, diagnoseSteps, statusName, STEP_MIN_FLUX } from './fba.js';
import { sampleFluxSpace, carbonCount } from './analysis_engine.js';
import { fluxMini, fluxNum } from './charts.js';

const esc = (s) => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const TOP_N = 10;               // pathways feasibility-tested automatically per run
const TOP_DEEP = 5;             // feasible pathways given automatic FVA + sampling
const SAMPLE_MIN = 10, SAMPLE_CAP = 200, SAMPLE_DEFAULT = 40;   // per-pathway browser caps, stated in the UI
const FLUX_EDGE_TOL = 1e-6;     // |flux| above this counts as carrying flux on the map

export async function initMode2(panel, ctx) {
  // Loading the engine and the media definitions decides whether Mode 2 exists
  // at all; a failure here disables the mode with the real reason.
  const [, mediaLib] = await Promise.all([getGLPK(), loadMedia()]);

  const state = {
    active: true,
    gemAcc: null, gem: null,
    customs: new Map(),          // label -> medium def (session only)
    currentLabel: null,
    working: null,               // {label, carbon_exchange, carbon_cap, components:{ex:lb}, edited}
    swapCarbon: true,
    sub: null, prod: null,
    run: null,                   // last feasibility run context
    runToken: 0,
    deep: new Map(),             // idx -> {fva, sample, pfba, rids} per-pathway FVA + sampling
    cancel: false,               // set by the Cancel control; checked between solves
    busyDeep: false,
  };

  // ---------- working-medium construction ----------
  function mediumDef(label) {
    return state.customs.get(label) || mediaLib[label] || null;
  }
  function loadWorking(label) {
    const def = mediumDef(label);
    if (!def) { state.working = null; state.currentLabel = null; return; }
    const components = {};
    for (const [ex, lb] of Object.entries(def.components || {})) {
      if (ex.startsWith('_')) continue;
      components[ex] = lb;
    }
    for (const [ex, lb] of Object.entries(def.supplements || {})) {
      if (ex.startsWith('_')) continue;
      if (!(ex in components)) components[ex] = lb;
    }
    state.working = {
      label,
      species: def.species || null,
      carbon_exchange: def.carbon_exchange || null,
      carbon_cap: def.carbon_cap ?? null,
      components,
      edited: false,
    };
    state.currentLabel = label;
  }

  function exToMet(exId) { return exId.replace(/^EX_/, ''); }
  function exName(exId) { return ctx.metName(exToMet(exId)) || ''; }

  // A GEM, medium or bound change makes every rendered feasibility claim
  // stale: clear the chips and details rather than letting them assert the
  // old setup.
  function invalidateRun(why) {
    if (!state.run) return;
    state.run = null;
    state.deep = new Map();
    state.runToken++;
    document.querySelectorAll('#results-body .feas-slot, #results-body .mode2-slot, #results-body .flux-slot')
      .forEach(el => { el.innerHTML = ''; });
    const s = document.querySelector('#results-summary .mode2-status');
    if (s) s.innerHTML = `<span class="status">${why || 'The GEM or medium changed'}; choose "Test the found pathways" in Simulate to re-test.</span>`;
    if (ctx.onSearchUpdate) ctx.onSearchUpdate('cleared');
  }

  // Bound edits and knockouts made in the GEM browser or the Analysis view
  // apply to every solve, so results rendered before them are stale.
  onEditsChanged((acc) => {
    if (acc === state.gemAcc) invalidateRun('Reaction bounds changed (edit or knockout)');
  });

  // Effective bounds for a run: the working components, with the substrate
  // swapped in as carbon source when enabled and the GEM has its exchange.
  function effectiveMedium() {
    const w = state.working;
    const out = { ...w.components };
    const info = { swapped: false, subEx: null, closedEx: null, note: '' };
    if (!state.swapCarbon || !state.sub) return { bounds: out, info };
    const subEx = substrateExchange();
    if (!subEx) {
      info.note = `The GEM has no exchange reaction for the substrate; the medium keeps ${w.carbon_exchange || 'its carbon source'} open.`;
      return { bounds: out, info };
    }
    if (w.carbon_exchange && w.carbon_exchange !== subEx) {
      out[w.carbon_exchange] = 0;
      info.closedEx = w.carbon_exchange;
    }
    out[subEx] = w.carbon_cap ?? -10;
    info.swapped = true;
    info.subEx = subEx;
    return { bounds: out, info };
  }

  function substrateExchange() {
    if (!state.gem || !state.sub) return null;
    const base = state.sub.replace(/_[a-z]+$/, '');
    const exId = `EX_${base}_e`;
    return state.gem.reactions.some(r => r.ex && r.id === exId) ? exId : null;
  }

  // ---------- panel ----------
  panel.innerHTML = `
    <div class="card mode2-card">
      <div class="mode2-head">
        <h2 id="m2-title">Flux context: GEM + medium</h2>
        <span class="status" id="m2-engine">GLPK 5.0 (WASM) loaded · solves in a background worker</span>
      </div>
      <details class="methodnote">
        <summary>Method</summary>
        <p>Feasibility is tested per pathway for one GEM on one medium: the LP maximises
        a demand on the product while every pathway step carries at least ${STEP_MIN_FLUX} mmol gDW<sup>-1</sup> h<sup>-1</sup>
        in the pathway direction. Exchanges not listed in the medium are closed. pFBA and FVA run
        per pathway on request. Definitions and limits: <a href="methods.html#analysis">methods page</a>.</p>
      </details>
      <div class="gem-toolbar">
        <div class="field">
          <label for="m2-gem">GEM (1 of 35)</label>
          <select id="m2-gem"><option value="">Choose a GEM…</option></select>
        </div>
        <div class="field">
          <label for="m2-medium">Medium (9 predefined)</label>
          <select id="m2-medium"><option value="">Choose a medium…</option></select>
        </div>
        <button class="btn small" id="m2-clone" type="button" disabled>Clone to custom medium</button>
        <span class="status" id="m2-status" role="status" aria-live="polite"></span>
      </div>
      <div id="m2-medium-editor" hidden></div>
      <div id="m2-swap" class="m2-swap" hidden></div>
    </div>`;

  const $p = (sel) => panel.querySelector(sel);
  const gemSel = $p('#m2-gem');
  const medSel = $p('#m2-medium');
  const statusEl = $p('#m2-status');

  // GEM options grouped by species
  {
    const bySp = new Map();
    for (const g of ctx.index.gems) {
      if (!bySp.has(g.species)) bySp.set(g.species, []);
      bySp.get(g.species).push(g);
    }
    for (const [sp, gems] of bySp) {
      const og = document.createElement('optgroup');
      og.label = `${sp} (${gems.length})`;
      for (const g of gems) {
        const o = document.createElement('option');
        o.value = g.acc;
        o.textContent = `${g.acc} · ${fmt.format(g.reactions)} rxns`;
        og.appendChild(o);
      }
      gemSel.appendChild(og);
    }
  }

  function rebuildMediumOptions() {
    const cur = medSel.value;
    medSel.innerHTML = '<option value="">Choose a medium…</option>';
    const og1 = document.createElement('optgroup');
    og1.label = `Predefined (${Object.keys(mediaLib).length})`;
    for (const label of Object.keys(mediaLib)) {
      const o = document.createElement('option');
      o.value = label; o.textContent = label;
      og1.appendChild(o);
    }
    medSel.appendChild(og1);
    if (state.customs.size) {
      const og2 = document.createElement('optgroup');
      og2.label = `Custom, this session (${state.customs.size})`;
      for (const label of state.customs.keys()) {
        const o = document.createElement('option');
        o.value = label; o.textContent = label;
        og2.appendChild(o);
      }
      medSel.appendChild(og2);
    }
    medSel.value = cur;
  }
  rebuildMediumOptions();

  gemSel.addEventListener('change', async () => {
    const acc = gemSel.value;
    invalidateRun();
    state.gem = null; state.gemAcc = null;
    setContext({ gem: acc || null });
    if (!acc) { renderMediumEditor(); renderSwapNote(); notifyState(); return; }
    try {
      statusEl.textContent = `Loading GEM ${acc} (about 1 MB)…`;
      const gem = await loadGem(acc);
      state.gem = gem; state.gemAcc = acc;
      statusEl.textContent = `${gem.species} · ${acc} · ${fmt.format(gem.reactions.length)} reactions`;
      ctx.setAccent(gem.species);
    } catch (e) {
      statusEl.textContent = `Could not load GEM ${acc}: ${e.message}. Choose it again to retry.`;
    }
    renderMediumEditor();
    renderSwapNote();
    notifyState();
  });

  medSel.addEventListener('change', () => {
    invalidateRun();
    loadWorking(medSel.value || null);
    setContext({ medium: medSel.value || null });
    $p('#m2-clone').disabled = !state.working;
    renderMediumEditor();
    renderSwapNote();
    notifyState();
  });

  function notifyState() { if (ctx.onStateChange) ctx.onStateChange(); }

  // Context written elsewhere (the Model stage select, the Engineer toolbar)
  // is applied here so the chosen GEM and medium carry into the flux setup.
  function applyContext(c) {
    if ((c.gem || '') !== gemSel.value) {
      const has = c.gem && ctx.index.gems.some(g => g.acc === c.gem);
      gemSel.value = has ? c.gem : '';
      gemSel.dispatchEvent(new Event('change'));
    }
    if ((c.medium || '') !== medSel.value) {
      const label = c.medium || '';
      const has = label && (mediaLib[label] || state.customs.has(label));
      medSel.value = has ? label : '';
      medSel.dispatchEvent(new Event('change'));
    }
  }
  onContext((c, changed) => {
    if (changed.includes('gem') || changed.includes('medium')) applyContext(c);
  });
  applyContext(getContext());

  $p('#m2-clone').addEventListener('click', () => {
    if (!state.working) return;
    const n = state.customs.size + 1;
    const label = `Custom ${n} · from ${state.working.label.replace(/^Custom \d+ · from /, '')}`;
    state.customs.set(label, {
      species: state.working.species,
      carbon_exchange: state.working.carbon_exchange,
      carbon_cap: state.working.carbon_cap,
      components: { ...state.working.components },
    });
    rebuildMediumOptions();
    medSel.value = label;
    loadWorking(label);
    setContext({ medium: label });
    renderMediumEditor();
    renderSwapNote();
    notifyState();
    statusEl.textContent = `Saved as "${label}" (held in browser memory for this session).`;
  });

  // ---------- medium editor ----------
  function renderMediumEditor() {
    const host = $p('#m2-medium-editor');
    const w = state.working;
    if (!w) { host.hidden = true; host.innerHTML = ''; return; }
    host.hidden = false;
    const exIds = Object.keys(w.components).sort((a, b) =>
      (a === w.carbon_exchange ? -1 : b === w.carbon_exchange ? 1 : a.localeCompare(b)));
    const gemExSet = state.gem ? new Set(state.gem.reactions.filter(r => r.ex).map(r => r.id)) : null;
    const unknown = gemExSet ? exIds.filter(e => !gemExSet.has(e)) : [];
    const spNote = (state.gem && w.species && state.gem.species && !state.gem.species.startsWith(w.species.split('_')[0]))
      ? `<p class="termination capped">Medium "${esc(w.label)}" is defined for ${esc(w.species.replace('_', ' '))}; the chosen GEM is ${esc(state.gem.species)}. Exchange ids may not match.</p>` : '';

    host.innerHTML = `
      <p class="count">${exIds.length} exchange components · carbon source
        ${w.carbon_exchange ? `<span class="mono">${esc(w.carbon_exchange)}</span> at lb ${w.carbon_cap ?? 'not recorded'}` : 'not recorded'}
        ${w.edited ? ' · edited (held in browser memory for this session)' : ''}
        ${unknown.length ? ` · ${unknown.length} of ${exIds.length} components have no exchange in ${esc(state.gemAcc)} and are ignored when solving` : ''}</p>
      ${spNote}
      <div class="tablewrap" style="max-width:640px;max-height:320px;overflow-y:auto">
        <table class="data" id="m2-medium-table">
          <thead><tr><th scope="col">Exchange</th><th scope="col">Metabolite</th>
            <th scope="col">lb (uptake, mmol gDW<sup>-1</sup> h<sup>-1</sup>)</th><th scope="col"></th></tr></thead>
          <tbody>${exIds.map(ex => `
            <tr data-ex="${esc(ex)}" class="${gemExSet && !gemExSet.has(ex) ? 'm2-unknown' : ''}">
              <td class="mono">${esc(ex)}${ex === w.carbon_exchange ? ' <span class="m2-tag">carbon source</span>' : ''}</td>
              <td>${esc(exName(ex))}</td>
              <td><input class="bound" type="number" step="any" value="${w.components[ex]}" aria-label="Lower bound of ${esc(ex)}"></td>
              <td><button class="btn small" type="button" data-remove="${esc(ex)}">Remove</button></td>
            </tr>`).join('')}
          </tbody>
        </table>
      </div>
      <div class="tablebar">
        <div class="field" style="flex:1 1 240px">
          <label for="m2-add-ex">Add exchange component${state.gem ? ` (from the ${fmt.format(state.gem.stats.exchanges)} exchanges of ${esc(state.gemAcc)})` : ''}</label>
          ${state.gem
            ? `<select id="m2-add-ex"><option value="">Choose an exchange…</option>${state.gem.reactions.filter(r => r.ex && !(r.id in w.components)).map(r => `<option value="${esc(r.id)}">${esc(r.id)}${exName(r.id) ? ' · ' + esc(exName(r.id)) : ''}</option>`).join('')}</select>`
            : '<select id="m2-add-ex" disabled><option>Choose a GEM to list its exchange reactions</option></select>'}
        </div>
        <button class="btn small" id="m2-add-btn" type="button" ${state.gem ? '' : 'disabled'}>Add at lb -1000</button>
      </div>`;

    host.querySelectorAll('input.bound').forEach(inp => {
      inp.addEventListener('change', () => {
        const ex = inp.closest('tr').dataset.ex;
        const v = parseFloat(inp.value);
        if (Number.isNaN(v)) return;
        invalidateRun();
        w.components[ex] = v;
        w.edited = true;
        if (ex === w.carbon_exchange) w.carbon_cap = v;
        renderMediumEditor();
        renderSwapNote();
      });
    });
    host.querySelectorAll('button[data-remove]').forEach(btn => {
      btn.addEventListener('click', () => {
        invalidateRun();
        delete w.components[btn.dataset.remove];
        if (btn.dataset.remove === w.carbon_exchange) { w.carbon_exchange = null; w.carbon_cap = null; }
        w.edited = true;
        renderMediumEditor();
        renderSwapNote();
      });
    });
    const addBtn = host.querySelector('#m2-add-btn');
    if (addBtn) addBtn.addEventListener('click', () => {
      const sel = host.querySelector('#m2-add-ex');
      if (!sel || !sel.value) return;
      invalidateRun();
      w.components[sel.value] = -1000;
      w.edited = true;
      renderMediumEditor();
      renderSwapNote();
    });
  }

  // ---------- substrate swap note ----------
  function renderSwapNote() {
    const host = $p('#m2-swap');
    const w = state.working;
    if (!w) { host.hidden = true; host.innerHTML = ''; return; }
    host.hidden = false;
    let line;
    if (!state.sub) {
      line = 'Pick a substrate above; when the search runs it replaces the medium\'s carbon source.';
    } else if (!state.gem) {
      line = 'Choose a GEM; the substrate swap needs the GEM\'s exchange reactions.';
    } else {
      const subEx = substrateExchange();
      const subNm = ctx.metName(state.sub) || state.sub;
      if (!state.swapCarbon) {
        line = `The substrate is not swapped in: the medium's carbon source ${w.carbon_exchange ? `<span class="mono">${esc(w.carbon_exchange)}</span>` : ''} stays open.`;
      } else if (!subEx) {
        line = `${esc(state.gemAcc)} has no exchange reaction for ${esc(subNm)} <span class="mono">${esc(state.sub)}</span>; the carbon source is not swapped and the medium runs as listed.`;
      } else if (subEx === w.carbon_exchange) {
        line = `The substrate ${esc(subNm)} already is the medium's carbon source (<span class="mono">${esc(subEx)}</span> at lb ${w.carbon_cap ?? 'not recorded'}).`;
      } else {
        line = `When the search runs, ${w.carbon_exchange ? `<span class="mono">${esc(w.carbon_exchange)}</span> closes and ` : ''}<span class="mono">${esc(subEx)}</span> (${esc(subNm)}) opens at lb ${w.carbon_cap ?? -10}, the medium's carbon cap.`;
      }
    }
    host.innerHTML = `
      <label class="mode" style="border-style:solid"><input type="checkbox" id="m2-swap-cb" ${state.swapCarbon ? 'checked' : ''}>
        Substrate replaces the medium's carbon source</label>
      <span class="status" aria-live="polite">${line}</span>`;
    host.querySelector('#m2-swap-cb').addEventListener('change', (e) => {
      invalidateRun();
      state.swapCarbon = e.target.checked;
      renderSwapNote();
    });
  }

  // ---------- feasibility run ----------
  function summaryLine() {
    let el = document.querySelector('#results-summary .mode2-status');
    if (!el) {
      el = document.createElement('div');
      el.className = 'mode2-status';
      ctx.resultsSummary().appendChild(el);
    }
    return el;
  }

  function chip(card, cls, text) {
    const slot = card && card.querySelector('.feas-slot');
    if (slot) slot.innerHTML = `<span class="feas-chip ${cls}">${text}</span>`;
  }

  function sampleCount() {
    const el = document.querySelector('#sim-samples');
    const v = el ? parseInt(el.value, 10) : SAMPLE_DEFAULT;
    if (Number.isNaN(v)) return SAMPLE_DEFAULT;
    return Math.min(Math.max(v, SAMPLE_MIN), SAMPLE_CAP);
  }

  function notifySearch(phase) { if (ctx.onSearchUpdate) ctx.onSearchUpdate(phase); }

  function feasibleIdxList() {
    if (!state.run) return [];
    return [...state.run.results.entries()]
      .filter(([, r]) => r.testable && r.feasible)
      .map(([i]) => i)
      .sort((a, b) => a - b);
  }

  async function runFeasibility(lastResults) {
    const token = ++state.runToken;
    state.cancel = false;
    const { res, prod } = lastResults;
    const pathways = res.pathways;
    state.run = null;
    state.deep = new Map();
    if (!pathways.length) return;
    const gem = state.gem, acc = state.gemAcc, w = state.working;
    const line = summaryLine();
    setCancelVisible(true);

    const target = productTarget(gem, prod);
    if (!target) {
      line.innerHTML = `<span class="status">${esc(acc)} does not contain the product <span class="mono">${esc(prod)}</span>; no pathway can be tested in this GEM.</span>`;
      setCancelVisible(false);
      return;
    }
    const { bounds, info } = effectiveMedium();
    state.run = { gem, acc, bounds, target, info, lastResults, results: new Map() };

    line.innerHTML = `<span class="status">Solving growth on ${esc(w.label)}…</span>`;
    let muText = 'not computed';
    try {
      const g = await maxGrowth(gem, acc, bounds);
      if (token !== state.runToken) return;
      muText = g.optimal ? `μmax ${g.mu.toFixed(3)} h<sup>-1</sup>` : `growth LP: ${esc(statusName(g.status))}`;
    } catch (e) {
      muText = `growth solve failed (${esc(e.message)})`;
    }
    const swapText = info.swapped
      ? `substrate in as carbon source (<span class="mono">${esc(info.subEx)}</span> at lb ${w.carbon_cap ?? -10}${info.closedEx ? `, <span class="mono">${esc(info.closedEx)}</span> closed` : ''})`
      : (info.note ? esc(info.note) : 'carbon source unchanged');
    state.run.muText = muText;
    state.run.swapText = swapText;

    const nTest = Math.min(TOP_N, pathways.length);
    let cancelledAt = null;
    for (let i = 0; i < nTest; i++) {
      if (token !== state.runToken) return;
      if (state.cancel) { cancelledAt = i; break; }
      line.innerHTML = `<span class="status">${esc(acc)} on ${esc(w.label)} · ${muText} · ${swapText} ·
        testing pathway feasibility ${i + 1} of ${nTest}${pathways.length > nTest ? ` (of ${fmt.format(pathways.length)} found; test the rest per card)` : ''}…</span>`;
      await testPathway(i, pathways[i], token);
      notifySearch('progress');
    }
    if (token !== state.runToken) return;

    // FVA + flux sampling for the shortest feasible pathways, deepest first
    const feasIdx = feasibleIdxList();
    const deepIdx = feasIdx.slice(0, TOP_DEEP);
    let deepDone = 0;
    for (const i of deepIdx) {
      if (token !== state.runToken) return;
      if (state.cancel) { if (cancelledAt == null) cancelledAt = nTest; break; }
      line.innerHTML = `<span class="status">${esc(acc)} on ${esc(w.label)} · FVA + flux sampling on feasible pathway
        #${i + 1} (${deepDone + 1} of ${deepIdx.length} deep runs; ${sampleCount()} samples each)…</span>`;
      await deepAnalyse(i, pathways[i], token);
      deepDone++;
      notifySearch('progress');
    }
    if (token !== state.runToken) return;
    setCancelVisible(false);

    const tested = state.run.results.size;
    const feas = feasibleIdxList().length;
    const cancelTxt = cancelledAt != null
      ? ` · cancelled by you after ${cancelledAt} of ${nTest} feasibility tests and ${deepDone} of ${deepIdx.length} deep runs; the remainder is untested, not infeasible` : '';
    line.innerHTML = `<span class="status">${esc(acc)} on ${esc(w.label)} · ${muText} · ${swapText} ·
      ${feas} of ${tested} tested pathways feasible${pathways.length > nTest ? ` (${fmt.format(pathways.length - nTest)} of ${fmt.format(pathways.length)} untested; use "Test this pathway" on a card)` : ''}
      · FVA + sampling on ${deepDone} of ${feas} feasible${cancelTxt}.</span>`;

    // already-rendered cards beyond the auto-tested set get an honest chip + button
    document.querySelectorAll('#results-body .pcard').forEach(c => {
      const i = +c.dataset.pwIdx;
      if (!state.run.results.has(i) && pathways[i]) {
        chip(c, 'na', 'not tested');
        renderCardDetail(i, pathways[i]);
      }
    });
    notifySearch('done');
  }

  // Testability comes from the GEM's own step realizations (fba.js), never
  // from the union presence mask: union ids can differ in stoichiometry per
  // GEM, so the mask can both over- and under-state carriage.
  async function testPathway(idx, pw, token) {
    const run = state.run;
    if (!run) return;
    const card = ctx.findCard(idx);
    chip(card, 'wait', 'solving…');
    try {
      const r = await pathwayFeasibility(run.gem, run.acc, run.bounds, pw, run.target);
      if (token !== state.runToken) return;
      run.results.set(idx, r);
      if (!r.testable) {
        chip(card, 'na', `not testable in ${esc(run.acc)} (no reaction realizes step ${r.missingStep + 1} of ${pw.len})`);
      } else if (r.feasible) {
        chip(card, 'ok', `feasible · ${r.productFlux.toFixed(2)} mmol gDW<sup>-1</sup> h<sup>-1</sup>`);
      } else {
        chip(card, 'bad', 'infeasible (this GEM + medium)');
      }
      renderCardDetail(idx, pw);
    } catch (e) {
      chip(card, 'na', `solve failed (${esc(e.message)})`);
    }
  }

  // ---------- yield vs substrate uptake, from the feasibility optimum ----------
  // Basis: the substrate exchange of the run (swapped-in substrate when active,
  // otherwise the medium's carbon source). mol/mol from the LP solution; C-mol
  // basis only when both formulas carry a carbon count. Rendered into the
  // pathway header's Yield fact.
  function yieldFactHTML(r) {
    const run = state.run;
    if (!r.vars || r.productFlux == null) return '<span class="fact-na">not computed</span>';
    const subEx = (run.info && run.info.subEx) || (state.working && state.working.carbon_exchange);
    if (!subEx) return '<span class="fact-na">not computed (no substrate exchange identified)</span>';
    const up = r.vars[subEx];
    if (up == null || up >= -1e-9) return `<span class="fact-na">not computed (no <span class="mono">${esc(subEx)}</span> uptake in the optimum)</span>`;
    const molmol = r.productFlux / Math.abs(up);
    const subMid = subEx.replace(/^EX_/, '');
    const prodMid = run.target.mid || run.lastResults.prod;
    const fS = run.gem.metabolites.find(m => m.id === subMid);
    const fP = run.gem.metabolites.find(m => m.id === prodMid) || run.gem.metabolites.find(m => m.id === run.lastResults.prod);
    const cS = carbonCount(fS && fS.formula), cP = carbonCount(fP && fP.formula);
    const cmol = (cS && cP) ? molmol * cP / cS : null;
    return `<strong>${molmol.toFixed(3)}</strong> mol/mol${cmol != null ? ` · ${cmol.toFixed(3)} C-mol/C-mol` : ''}
      <span class="fact-note">${cmol == null ? 'C-mol basis not computed: formula lacks a carbon count · ' : ''}basis: product flux vs <span class="mono">${esc(subEx)}</span> uptake ${Math.abs(up).toFixed(3)} mmol gDW<sup>-1</sup> h<sup>-1</sup></span>`;
  }

  // A header fact of the expanded panel (Feasibility, Max product flux, Yield).
  function setFact(card, key, html) {
    const dd = card && card.querySelector(`.pwh-facts [data-pwh="${key}"] dd`);
    if (dd) dd.innerHTML = html;
  }
  const FACT_NA = '<span class="fact-na">not computed</span>';

  // ---------- FVA + flux sampling per pathway (the deep run) ----------
  // FVA: min and max flux of each pathway reaction with the product held at
  // >= 99% of the pathway optimum and every step constraint active (fba.js).
  // Sampling: random-objective vertex sampling of the SAME constrained space
  // (analysis_engine), so the sampled distribution sits inside the FVA
  // envelope and shows where within it flux typically falls.
  const FVA_RID_CAP = 30;

  function pathwayRids(pw) {
    const rids = [];
    for (const vars of stepRealizations(state.run.gem, pw)) {
      for (const v of vars) if (!rids.includes(v.name)) rids.push(v.name);
    }
    return rids;
  }

  async function samplePathway(pw, productOpt, n, onProgress) {
    const run = state.run;
    const glpk = await getGLPK();
    const { stepCons } = stepConstraints(run.gem, pw);
    if (!stepCons) return { ok: false, reason: 'pathway not carried by this GEM' };
    const lp = buildLP(glpk, run.gem, run.bounds, {
      acc: run.acc, extraCols: run.target.extraCols, stepCons,
      objective: { direction: 'max', vars: [] },
    });
    for (const b of lp.bounds) if (b.name === run.target.id) {
      b.lb = productOpt * 0.99;
      b.ub = Math.max(b.ub, productOpt);
      b.type = boundType(glpk, b.lb, b.ub);
    }
    const S = {
      glpk, gem: run.gem, acc: run.acc, mediumBounds: run.bounds, lp,
      bIdx: new Map(lp.bounds.map(b => [b.name, b])),
      rowIdx: new Map(lp.subjectTo.map(row => [row.name, row])),
      target: run.target, biomass: run.gem.stats.biomass_id,
      rxnById: new Map(run.gem.reactions.map(r => [r.id, r])),
      kos: new Set(), solves: 0,
    };
    return sampleFluxSpace(S, { n, biomassFrac: 0, seed: 7, onProgress, shouldStop: () => state.cancel });
  }

  async function deepAnalyse(idx, pw, token) {
    const run = state.run;
    if (!run) return;
    const r = run.results.get(idx);
    if (!r || !r.testable || !r.feasible) return;
    const card = ctx.findCard(idx);
    const prog = card && card.querySelector('[data-m2prog]');
    const say = (t) => { if (prog) prog.textContent = t; };
    state.busyDeep = true;
    try {
      // pFBA: one parsimonious distribution at the product optimum
      say('Solving pFBA…');
      const p = await pathwayPFBA(run.gem, run.acc, run.bounds, pw, run.target, r.productFlux);
      if (token !== state.runToken) return;
      const allRids = pathwayRids(pw);
      const rids = allRids.slice(0, FVA_RID_CAP);
      // FVA over the pathway reactions
      const fva = await pathwayFVA(run.gem, run.acc, run.bounds, pw, run.target, r.productFlux, rids,
        (d, t) => say(`FVA ${d} of ${t} pathway reactions…`), () => state.cancel);
      if (token !== state.runToken) return;
      // flux sampling of the same constrained space
      const n = sampleCount();
      say(`Flux sampling 0 of ${n}…`);
      const sample = await samplePathway(pw, r.productFlux, n,
        (i, t) => say(`Flux sampling ${i} of ${t}…`));
      if (token !== state.runToken) return;
      state.deep.set(idx, {
        rids, allRids, fva, sample,
        pfba: p.optimal ? p.fluxes : null, pfbaStatus: p.status,
        stamp: { acc: run.acc, medium: state.working.label, samples: sample.ok ? sample.samples : 0, requested: n },
      });
      fillFluxSlots(idx, pw);
      renderCardDetail(idx, pw);
      say('');
    } catch (e) {
      say(`FVA + sampling failed (${e.message}).`);
    } finally {
      state.busyDeep = false;
    }
  }

  // Shared view over a deep run: one flux domain for every mini chart of the
  // pathway (so the bars compare across steps) and the reaction with the
  // narrowest FVA envelope. Computed once per deep result, from stored solves.
  function deepView(d) {
    if (d.view) return d.view;
    const vals = [];
    let limRid = null, limW = Infinity, nRanged = 0;
    for (const rid of d.rids) {
      const rr = d.fva && d.fva.ranges[rid];
      if (rr && rr.min != null && rr.max != null) {
        vals.push(rr.min, rr.max);
        nRanged++;
        const w = rr.max - rr.min;
        if (w < limW) { limW = w; limRid = rid; }
      }
      const ss = d.sample && d.sample.ok && d.sample.stats.get(rid);
      if (ss) vals.push(ss.min, ss.max);
      const v = d.pfba ? d.pfba[rid] : null;
      if (v != null) vals.push(v);
    }
    if (!vals.length) { d.view = { lo: null, hi: null, limRid: null, limW: null, nRanged: 0 }; return d.view; }
    let lo = Math.min(...vals, 0), hi = Math.max(...vals, 1e-9);
    const pad = (hi - lo) * 0.05 || 0.5;
    d.view = { lo: lo < 0 ? lo - pad : lo, hi: hi + pad, limRid, limW, nRanged };
    return d.view;
  }

  // Per-reaction flux panel in the step rows: an FVA envelope + sampled
  // distribution mini on the shared domain, the numbers with units, and a tag
  // on the reaction with the narrowest FVA range.
  function fillFluxSlots(idx, pw) {
    const run = state.run;
    const card = ctx.findCard(idx);
    const d = state.deep.get(idx);
    if (!run || !card || !d) return;
    const view = deepView(d);
    const gemRxns = new Set(run.gem.reactions.map(x => x.id));
    card.querySelectorAll('.rxn-alt').forEach(row => {
      const rid = row.dataset.rxn;
      const fslot = row.querySelector('.flux-slot');
      if (!fslot) return;
      if (!gemRxns.has(rid)) { fslot.innerHTML = `<span class="fluxna">not in ${esc(run.acc)}</span>`; return; }
      if (!d.rids.includes(rid)) {
        fslot.innerHTML = `<span class="fluxna">not computed: beyond the FVA cap (${d.rids.length} of ${d.allRids.length} pathway reactions ranged)</span>`;
        return;
      }
      const rr = d.fva && d.fva.ranges[rid];
      const fva = (rr && rr.min != null && rr.max != null) ? rr : null;
      const ss = (d.sample && d.sample.ok && d.sample.stats.get(rid)) || null;
      const v = d.pfba ? (d.pfba[rid] ?? null) : null;
      const lim = (view.limRid === rid && view.nRanged > 1)
        ? '<span class="limit-tag" title="Smallest max minus min FVA flux among this pathway&#39;s ranged reactions">narrowest FVA range</span>' : '';
      const mini = view.lo != null ? fluxMini(fva, ss, v, { lo: view.lo, hi: view.hi, width: 240 }) : '';
      fslot.innerHTML = `${lim}${mini}
        <div class="fluxnums">
          <span class="fn-fva">FVA ${fva ? `[${fluxNum(fva.min)}, ${fluxNum(fva.max)}]` : 'not computed'}</span>
          <span class="fn-p">pFBA ${v != null ? fluxNum(v) : 'not computed'}</span>
          <span class="fn-s">${ss ? `med ${fluxNum(ss.median)} · 5-95% [${fluxNum(ss.p5)}, ${fluxNum(ss.p95)}]` : 'sampling not computed'}</span>
          <span class="fn-u">mmol gDW<sup>-1</sup> h<sup>-1</sup></span>
        </div>`;
    });
  }

  // Per-card Simulate details. The verdict, max product flux and yield go to
  // the pathway header facts; the panel's Simulate strip carries the method
  // caption with its denominators, the chart legend, and the run controls.
  function renderCardDetail(idx, pw) {
    const run = state.run;
    const card = ctx.findCard(idx);
    if (!run || !card) return;
    const slot = card.querySelector('.mode2-slot');
    if (!slot) return;
    const medium = state.working ? state.working.label : '';
    const onLine = `<span class="fact-note">on <span class="mono">${esc(run.acc)}</span> · ${esc(medium)}</span>`;
    const r = run.results.get(idx);
    if (!r) {
      setFact(card, 'feas', `<span class="fact-na">not tested on <span class="mono">${esc(run.acc)}</span> · ${esc(medium)}</span>`);
      setFact(card, 'flux', FACT_NA);
      setFact(card, 'yield', FACT_NA);
      slot.innerHTML = `<div class="m2-detail"><span class="status">Feasibility not tested on this GEM + medium.</span>
        <button class="btn small" type="button" data-m2test>Test this pathway</button>
        <span class="status" data-m2prog role="status" aria-live="polite"></span></div>`;
      slot.querySelector('[data-m2test]').addEventListener('click', async () => {
        await testPathway(idx, pw, state.runToken);
        notifySearch('done');
      });
      return;
    }
    if (!r.testable) {
      setFact(card, 'feas', `<span class="feas-chip na">not testable in ${esc(run.acc)}</span>`);
      setFact(card, 'flux', FACT_NA);
      setFact(card, 'yield', FACT_NA);
      slot.innerHTML = `<div class="m2-detail"><span class="status">Not testable: no reaction of ${esc(run.acc)} interconverts the metabolite pair of step ${r.missingStep != null ? r.missingStep + 1 : '?'} of ${pw.len}.</span></div>`;
      return;
    }
    if (!r.feasible) {
      setFact(card, 'feas', `<span class="feas-chip bad">infeasible</span> ${onLine}`);
      setFact(card, 'flux', `<span class="fact-na">none: LP ${esc(statusName(r.status))}</span>`);
      setFact(card, 'yield', FACT_NA);
      slot.innerHTML = `<div class="m2-detail">
        <span class="status">Infeasible: under ${esc(run.acc)}'s bounds on ${esc(medium)}, no steady-state flux carries every step at ≥ ${STEP_MIN_FLUX} while producing the product (LP: ${esc(statusName(r.status))}).</span>
        <button class="btn small" type="button" data-m2diag>Diagnose steps</button>
        <span class="status" data-m2prog role="status" aria-live="polite"></span>
      </div>`;
      slot.querySelector('[data-m2diag]').addEventListener('click', () => diagnose(idx, pw));
      return;
    }

    setFact(card, 'feas', `<span class="feas-chip ok">feasible</span> ${onLine}`);
    setFact(card, 'flux', `<strong>${r.productFlux.toFixed(3)}</strong> mmol gDW<sup>-1</sup> h<sup>-1</sup>
      <span class="fact-note">product export maximum; every step at ≥ ${STEP_MIN_FLUX} mmol gDW<sup>-1</sup> h<sup>-1</sup></span>`);
    setFact(card, 'yield', yieldFactHTML(r));

    const d = state.deep.get(idx);
    if (!d) {
      slot.innerHTML = `<div class="m2-detail">
        <span class="status">Per-reaction flux ranges not computed yet for this pathway.</span>
        <button class="btn small" type="button" data-m2deep>Run FVA + flux sampling (${sampleCount()} samples)</button>
        <span class="status" data-m2prog role="status" aria-live="polite"></span></div>`;
      slot.querySelector('[data-m2deep]').addEventListener('click', async () => {
        state.cancel = false;
        await deepAnalyse(idx, pw, state.runToken);
        notifySearch('done');
      });
      return;
    }

    const view = deepView(d);
    const sNote = d.sample && d.sample.ok
      ? `${d.sample.samples} of ${d.stamp.requested} requested samples solved${d.sample.failed ? ` (${d.sample.failed} failed)` : ''}${d.sample.cancelled ? '; sampling cancelled early' : ''}, seed ${d.sample.seed}`
      : `sampling not computed${d.sample && d.sample.reason ? ` (${esc(d.sample.reason)})` : ''}`;
    const fvaNote = d.fva && d.fva.cancelled
      ? `FVA cancelled after ${d.fva.done} of ${d.fva.total} reactions` :
      `FVA over ${d.rids.length}${d.allRids.length > d.rids.length ? ` of ${d.allRids.length}` : ''} pathway reaction${d.allRids.length > 1 ? 's' : ''} in ${esc(run.acc)}`;
    const limNote = (view.limRid && view.nRanged > 1)
      ? `<span class="status">Narrowest FVA range of the ${view.nRanged} reactions ranged: <span class="mono">${esc(view.limRid)}</span> (width ${fluxNum(view.limW)} mmol gDW<sup>-1</sup> h<sup>-1</sup>); tagged in its step.</span>` : '';
    slot.innerHTML = `<div class="m2-deepinfo">
        <span class="status">Flux panels beside each step: ${fvaNote}, product held at ≥ 99% of this pathway's optimum with every step constraint active.
          Sampling: random-objective vertex sampling of the same constrained space; ${sNote}.
          All panels share one flux axis, ${fluxNum(view.lo)} to ${fluxNum(view.hi)} mmol gDW<sup>-1</sup> h<sup>-1</sup>.</span>
        ${limNote}
        <div class="m2-legend">
          <span class="lgi"><span class="lg-env" aria-hidden="true"></span>FVA range</span>
          <span class="lgi"><span class="lg-pfba" aria-hidden="true"></span>pFBA optimum</span>
          <span class="lgi"><span class="lg-band" aria-hidden="true"></span>sampled 5-95%</span>
          <span class="lgi"><span class="lg-med" aria-hidden="true"></span>sampled median</span>
          <span class="lgi">whisker = sampled min to max</span>
        </div>
      </div>
      <div class="m2-detail">
        <button class="btn small" type="button" data-m2deep>Re-run FVA + sampling (${sampleCount()} samples)</button>
        <span class="status" data-m2prog role="status" aria-live="polite"></span>
      </div>`;
    slot.querySelector('[data-m2deep]').addEventListener('click', async () => {
      state.cancel = false;
      await deepAnalyse(idx, pw, state.runToken);
      notifySearch('done');
    });
  }

  // Per-step relaxation: which steps fail alone, measured, never guessed.
  async function diagnose(idx, pw) {
    const run = state.run;
    const card = ctx.findCard(idx);
    if (!run || !card) return;
    const prog = card.querySelector('[data-m2prog]');
    const btn = card.querySelector('[data-m2diag]');
    if (btn) btn.disabled = true;
    const token = state.runToken;
    try {
      prog.textContent = `Solving ${pw.len} single-step relaxations…`;
      const diag = await diagnoseSteps(run.gem, run.acc, run.bounds, pw, run.target);
      if (token !== state.runToken || !diag) return;
      const bad = diag.filter(d => !d.ok);
      if (bad.length) {
        prog.innerHTML = bad.map(d => {
          const st = pw.steps[d.step];
          return `Step ${d.step + 1} (<span class="mono">${esc(st.from)}</span> → <span class="mono">${esc(st.to)}</span>) cannot carry flux in the pathway direction under this GEM's bounds and medium.`;
        }).join('<br>') + `<br>${bad.length} of ${pw.len} steps ${bad.length === 1 ? 'fails' : 'fail'} alone.`;
      } else {
        prog.textContent = `Each of the ${pw.len} steps is feasible alone; the combination of all ${pw.len} step constraints is what fails.`;
      }
    } catch (e) {
      prog.textContent = `Diagnosis failed (${e.message}).`;
    } finally {
      if (btn) btn.disabled = false;
    }
  }

  function setCancelVisible(v) {
    const btn = document.querySelector('#sim-cancel');
    if (btn) btn.hidden = !v;
  }

  // Re-decorate any pathway card the list renders or re-renders ("Show more"
  // batches, filter and sort re-renders): tested cards get their stored chip
  // and detail back, untested ones an honest "not tested" + test button.
  function decorateCard(node) {
    if (!state.run) return;
    const idx = +node.dataset.pwIdx;
    const pw = state.run.lastResults.res.pathways[idx];
    if (!pw) return;
    const r = state.run.results.get(idx);
    if (!r) {
      chip(node, 'na', 'not tested');
    } else if (!r.testable) {
      chip(node, 'na', `not testable in ${esc(state.run.acc)}`);
    } else if (r.feasible) {
      chip(node, 'ok', `feasible · ${r.productFlux.toFixed(2)} mmol gDW<sup>-1</sup> h<sup>-1</sup>`);
    } else {
      chip(node, 'bad', 'infeasible (this GEM + medium)');
    }
    renderCardDetail(idx, pw);
    if (state.deep.has(idx)) fillFluxSlots(idx, pw);
  }

  const mo = new MutationObserver((muts) => {
    if (!state.run || !state.active) return;
    for (const m of muts) for (const n of m.addedNodes) {
      if (n.nodeType === 1 && n.classList && n.classList.contains('pcard')) decorateCard(n);
    }
  });
  const resultsBody = document.querySelector('#results-body');
  if (resultsBody) mo.observe(resultsBody, { childList: true, subtree: true });

  // ---------- public api ----------
  return {
    setActive(v) { state.active = v; },
    notReadyReason() {
      if (!state.gemAcc || !state.gem) return 'Choose a GEM (1 of 35) above; the choice carries across every stage.';
      if (!state.working) return 'Choose a medium (9 predefined, or a custom clone) above.';
      return null;
    },
    onEndpointsChange(sub, prod) {
      state.sub = sub; state.prod = prod;
      renderSwapNote();
    },
    runFeasibility,
    cancelRun() { state.cancel = true; },
    // Read-only view of the current run for the results list (sorting,
    // filtering, the length histogram and the map tiers).
    getRun() {
      if (!state.run) return null;
      const feasibleIdx = feasibleIdxList();
      const bestIdx = feasibleIdx.length ? feasibleIdx[0] : null;
      let fluxRids = null, bestWeights = null;
      if (bestIdx != null && state.deep.has(bestIdx)) {
        const d = state.deep.get(bestIdx);
        if (d.pfba) {
          fluxRids = Object.entries(d.pfba)
            .filter(([, v]) => v != null && Math.abs(v) > FLUX_EDGE_TOL)
            .map(([rid]) => rid);
          const pw = state.run.lastResults.res.pathways[bestIdx];
          bestWeights = stepRealizations(state.run.gem, pw).map(vars => {
            let s = 0;
            for (const v of vars) {
              const f = d.pfba[v.name];
              if (f == null) continue;
              s += v.coef * f;
            }
            return Math.abs(s);
          });
        }
      }
      return {
        acc: state.run.acc,
        medium: state.working ? state.working.label : null,
        results: state.run.results,
        deep: state.deep,
        feasibleIdx, bestIdx, fluxRids, bestWeights,
        total: state.run.lastResults.res.pathways.length,
      };
    },
    gemObject() { return state.gem; },
  };
}
