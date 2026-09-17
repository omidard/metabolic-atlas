// Mode 2: flux feasibility of Mode-1 pathways for a chosen GEM on an editable
// medium. Owns the media panel, the substrate-as-carbon-source swap, and the
// per-pathway feasibility / pFBA / FVA runs. All solving happens in the GLPK
// worker; every result carries its denominator; absent stays absent.

import { loadGem, loadMedia, fmt } from './data.js';
import { getGLPK, maxGrowth, productTarget, pathwayFeasibility, pathwayPFBA, pathwayFVA, stepConstraints, diagnoseSteps, statusName, STEP_MIN_FLUX } from './fba.js';

const esc = (s) => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const TOP_N = 10;               // pathways tested automatically per run

export async function initMode2(panel, ctx) {
  // Loading the engine and the media definitions decides whether Mode 2 exists
  // at all; a failure here disables the mode with the real reason.
  const [, mediaLib] = await Promise.all([getGLPK(), loadMedia()]);

  const state = {
    active: false,
    gemAcc: null, gem: null,
    customs: new Map(),          // label -> medium def (session only)
    currentLabel: null,
    working: null,               // {label, carbon_exchange, carbon_cap, components:{ex:lb}, edited}
    swapCarbon: true,
    sub: null, prod: null,
    run: null,                   // last feasibility run context
    runToken: 0,
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

  // A GEM or medium change makes every rendered feasibility claim stale:
  // clear the chips and details rather than letting them assert the old setup.
  function invalidateRun() {
    if (!state.run) return;
    state.run = null;
    state.runToken++;
    document.querySelectorAll('#results-body .feas-slot, #results-body .mode2-slot, #results-body .flux-slot')
      .forEach(el => { el.innerHTML = ''; });
    const s = document.querySelector('#results-summary .mode2-status');
    if (s) s.innerHTML = '<span class="status">The GEM or medium changed; run the search again to re-test feasibility.</span>';
  }

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
        <h2 id="m2-title">Mode 2 · flux feasibility</h2>
        <span class="status" id="m2-engine">GLPK 5.0 (WASM) loaded · solves in a background worker</span>
      </div>
      <p class="sub">Feasibility is tested per pathway for one GEM on one medium: the LP maximises
      a demand on the product while every pathway step carries at least ${STEP_MIN_FLUX} mmol gDW<sup>-1</sup> h<sup>-1</sup>
      in the pathway direction. Exchanges not listed in the medium are closed. pFBA and FVA run per pathway on request.</p>
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
    if (!acc) { renderMediumEditor(); renderSwapNote(); return; }
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
  });

  medSel.addEventListener('change', () => {
    invalidateRun();
    loadWorking(medSel.value || null);
    $p('#m2-clone').disabled = !state.working;
    renderMediumEditor();
    renderSwapNote();
  });

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
    renderMediumEditor();
    renderSwapNote();
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

  async function runFeasibility(lastResults) {
    const token = ++state.runToken;
    const { res, prod } = lastResults;
    const pathways = res.pathways;
    state.run = null;
    if (!pathways.length) return;
    const gem = state.gem, acc = state.gemAcc, w = state.working;
    const line = summaryLine();

    const target = productTarget(gem, prod);
    if (!target) {
      line.innerHTML = `<span class="status">${esc(acc)} does not contain the product <span class="mono">${esc(prod)}</span>; no pathway can be tested in this GEM.</span>`;
      return;
    }
    const { bounds, info } = effectiveMedium();
    state.run = { gem, acc, bounds, target, lastResults, results: new Map() };

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

    const nTest = Math.min(TOP_N, pathways.length);
    for (let i = 0; i < nTest; i++) {
      if (token !== state.runToken) return;
      line.innerHTML = `<span class="status">${esc(acc)} on ${esc(w.label)} · ${muText} · ${swapText} ·
        testing pathway feasibility ${i + 1} of ${nTest}${pathways.length > nTest ? ` (of ${fmt.format(pathways.length)} pathways; test the rest per card)` : ''}…</span>`;
      await testPathway(i, pathways[i], token);
    }
    if (token !== state.runToken) return;
    const feas = [...state.run.results.values()].filter(r => r.testable && r.feasible).length;
    line.innerHTML = `<span class="status">${esc(acc)} on ${esc(w.label)} · ${muText} · ${swapText} ·
      ${feas} of ${nTest} tested pathways feasible${pathways.length > nTest ? ` (${fmt.format(pathways.length - nTest)} more untested; use "Test feasibility" on a card)` : ''}.</span>`;

    // already-rendered cards beyond the auto-tested set get an honest chip + button
    document.querySelectorAll('#results-body .pcard').forEach(c => {
      const i = +c.dataset.pwIdx;
      if (!state.run.results.has(i) && pathways[i]) {
        chip(c, 'na', 'not tested');
        renderCardDetail(i, pathways[i]);
      }
    });
  }

  async function testPathway(idx, pw, token) {
    const run = state.run;
    if (!run) return;
    const card = ctx.findCard(idx);
    const accIdx = ctx.graphMeta.accs.findIndex(a => a.acc === run.acc);
    const carried = accIdx >= 0 && ((pw.mask >> BigInt(accIdx)) & 1n) === 1n;
    if (!carried) {
      const { missingStep } = stepConstraints(run.gem, pw);
      const stepTxt = missingStep != null ? ` (missing step ${missingStep + 1} of ${pw.len})` : '';
      run.results.set(idx, { testable: false });
      chip(card, 'na', `not carried by ${esc(run.acc)}${stepTxt}`);
      renderCardDetail(idx, pw);
      return;
    }
    chip(card, 'wait', 'solving…');
    try {
      const r = await pathwayFeasibility(run.gem, run.acc, run.bounds, pw, run.target);
      if (token !== state.runToken) return;
      run.results.set(idx, r);
      if (!r.testable) {
        chip(card, 'na', `not carried by ${esc(run.acc)} (missing step ${r.missingStep + 1} of ${pw.len})`);
      } else if (r.feasible) {
        chip(card, 'ok', `feasible · product ${r.productFlux.toFixed(2)} mmol gDW<sup>-1</sup> h<sup>-1</sup>`);
      } else {
        chip(card, 'bad', 'infeasible (this GEM + medium)');
      }
      renderCardDetail(idx, pw);
    } catch (e) {
      chip(card, 'na', `solve failed (${esc(e.message)})`);
    }
  }

  // Per-card Mode-2 details: definition line + pFBA/FVA button + flux fill-in.
  function renderCardDetail(idx, pw) {
    const run = state.run;
    const card = ctx.findCard(idx);
    if (!run || !card) return;
    const slot = card.querySelector('.mode2-slot');
    if (!slot) return;
    const r = run.results.get(idx);
    if (!r) {
      slot.innerHTML = `<div class="m2-detail"><span class="status">Not tested on this medium.</span>
        <button class="btn small" type="button" data-m2test>Test feasibility</button></div>`;
      slot.querySelector('[data-m2test]').addEventListener('click', () => testPathway(idx, pw, state.runToken));
      return;
    }
    if (!r.testable) {
      slot.innerHTML = `<div class="m2-detail"><span class="status">Not testable: ${esc(run.acc)} does not carry every step of this pathway.</span></div>`;
      return;
    }
    if (!r.feasible) {
      slot.innerHTML = `<div class="m2-detail">
        <span class="status">Infeasible: under ${esc(run.acc)}'s bounds on ${esc(state.working.label)}, no steady-state flux carries every step at ≥ ${STEP_MIN_FLUX} while producing the product (LP: ${esc(statusName(r.status))}).</span>
        <button class="btn small" type="button" data-m2diag>Diagnose steps</button>
        <span class="status" data-m2prog role="status" aria-live="polite"></span>
      </div>`;
      slot.querySelector('[data-m2diag]').addEventListener('click', () => diagnose(idx, pw));
      return;
    }
    slot.innerHTML = `<div class="m2-detail">
      <span class="status">Feasible on ${esc(state.working.label)}: max product flux ${r.productFlux.toFixed(3)} mmol gDW<sup>-1</sup> h<sup>-1</sup> with every step at ≥ ${STEP_MIN_FLUX}.</span>
      <button class="btn small" type="button" data-m2flux>pFBA + FVA on this pathway</button>
      <span class="status" data-m2prog role="status" aria-live="polite"></span>
    </div>`;
    slot.querySelector('[data-m2flux]').addEventListener('click', () => fluxDetail(idx, pw, r));
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

  async function fluxDetail(idx, pw, feas) {
    const run = state.run;
    const card = ctx.findCard(idx);
    if (!run || !card) return;
    const prog = card.querySelector('[data-m2prog]');
    const btn = card.querySelector('[data-m2flux]');
    if (btn) btn.disabled = true;
    const token = state.runToken;
    try {
      prog.textContent = 'Solving pFBA…';
      const p = await pathwayPFBA(run.gem, run.acc, run.bounds, pw, run.target, feas.productFlux);
      if (token !== state.runToken) return;
      if (!p.optimal) { prog.textContent = `pFBA: ${statusName(p.status)}; no flux distribution to show.`; return; }

      // reactions to report: the GEM-carried alternatives of each step
      const gemRxns = new Set(run.gem.reactions.map(x => x.id));
      const rids = [];
      for (const st of pw.steps) for (const alt of st.rxns) {
        if (gemRxns.has(alt.id) && !rids.includes(alt.id)) rids.push(alt.id);
      }
      const fva = await pathwayFVA(run.gem, run.acc, run.bounds, pw, run.target, feas.productFlux, rids,
        (d, t) => { prog.textContent = `FVA ${d} of ${t} reactions…`; });
      if (token !== state.runToken) return;

      // fill per-reaction flux + FVA into the step rows
      card.querySelectorAll('.rxn-alt').forEach(row => {
        const rid = row.dataset.rxn;
        const fslot = row.querySelector('.flux-slot');
        if (!fslot) return;
        if (!gemRxns.has(rid)) { fslot.innerHTML = `<span class="fluxval na">not in ${esc(run.acc)}</span>`; return; }
        const v = p.fluxes[rid];
        const rr = fva.ranges[rid];
        const vTxt = v == null ? 'not computed' : v.toFixed(3);
        const rTxt = rr && rr.min != null && rr.max != null ? ` · FVA [${rr.min.toFixed(3)}, ${rr.max.toFixed(3)}]` : '';
        fslot.innerHTML = `<span class="fluxval">v ${vTxt}${rTxt}</span>`;
      });

      // per-step directed flux magnitude drives edge thickness on the map
      const weights = pw.steps.map(st => {
        let s = 0;
        for (const alt of st.rxns) {
          const v = p.fluxes[alt.id];
          if (v == null) continue;
          s += (alt.dir === 'rev' ? -1 : 1) * v;
        }
        return Math.abs(s);
      });
      ctx.showFluxOnMap(pw, idx, weights);
      prog.innerHTML = `pFBA at product ${feas.productFlux.toFixed(3)} mmol gDW<sup>-1</sup> h<sup>-1</sup>; FVA over ${rids.length} pathway reactions at ≥ 99% of the product optimum. Edge thickness on the map scales with each step's |flux|.`;
    } catch (e) {
      prog.textContent = `Flux solve failed (${e.message}).`;
    } finally {
      if (btn) btn.disabled = false;
    }
  }

  // Decorate pathway cards rendered after a run ("Show more" batches).
  const mo = new MutationObserver((muts) => {
    if (!state.run || !state.active) return;
    for (const m of muts) for (const n of m.addedNodes) {
      if (n.nodeType === 1 && n.classList && n.classList.contains('pcard')) {
        const idx = +n.dataset.pwIdx;
        const pw = state.run.lastResults.res.pathways[idx];
        if (pw && !state.run.results.has(idx)) {
          chip(n, 'na', 'not tested');
          renderCardDetail(idx, pw);
        }
      }
    }
  });
  const resultsBody = document.querySelector('#results-body');
  if (resultsBody) mo.observe(resultsBody, { childList: true, subtree: true });

  // ---------- public api ----------
  return {
    setActive(v) { state.active = v; panel.hidden = !v; },
    notReadyReason() {
      if (!state.gemAcc || !state.gem) return 'Mode 2 needs a GEM: choose one of the 35 models in the Mode 2 panel.';
      if (!state.working) return 'Mode 2 needs a medium: choose one of the 9 predefined media (or a custom clone) in the Mode 2 panel.';
      return null;
    },
    onEndpointsChange(sub, prod) {
      state.sub = sub; state.prod = prod;
      renderSwapNote();
    },
    runFeasibility,
  };
}
