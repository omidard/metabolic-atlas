// Analysis view: single-reaction knockouts, shadow-price and reduced-cost
// estimates, greedy growth-coupling search, production envelope. All numbers
// come from GLPK solves in the browser on the chosen GEM + medium, respecting
// GEM-browser bound edits. Absent renders as "not computed"; every count
// carries its denominator; a heuristic that fails says so and never claims
// impossibility.

import { loadGem, loadMedia, setEdit, clearEdit, listEdits, onEditsChanged, fmt, downloadBlob, csvEscape } from './data.js';
import { statusName } from './fba.js';
import {
  analysisTarget, makeSession, baseState, koSweep, sweepScope,
  shadowPrices, reducedCosts, couplingSearch, productionEnvelope, floorYield,
  ZERO_MU, FLUX_TOL, COSTLY_FRAC, SHADOW_EPS,
} from './analysis_engine.js';

const esc = (s) => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
// Numbers within 1e-9 of zero print as zero (solver tolerance; methods page).
const fnum = (x, d = 4) => (x == null ? 'not computed' : (Math.abs(x) < 1e-9 ? (0).toFixed(d) : x.toFixed(d)));
const UNIT = 'mmol gDW<sup>-1</sup> h<sup>-1</sup>';

const CLS_LABEL = {
  lethal: 'lethal', costly: 'costly', neutral: 'neutral',
  beneficial: 'raises product floor', unsolved: 'unsolved',
};
const SP_BATCH = 200;    // shadow-price metabolites per computed batch
const RC_BATCH = 100;    // reduced-cost reactions per computed batch

export async function initAnalysis(root, ctx) {
  const mediaLib = await loadMedia();

  const st = {
    acc: null, gem: null,
    mediumLabel: null,
    sub: null, prod: null, swapCarbon: true,
    frac: 0.9,
    base: null, baseToken: 0,
    busy: null,              // name of the running long analysis, or null
    cancelFlag: false,
    editStamp: 0,            // bumped on every bound change
    sweep: null,             // {rows, refs, stamp, cancelled, sweptOf, ...}
    search: null,            // coupling result {res, stamp}
    env: null,               // {wt:{...}, ko:{...}|null, stamp}
    sp: new Map(),           // computed shadow prices mid -> rec (display cache)
    rc: new Map(),
    spStampNote: null,
  };

  onEditsChanged((acc) => {
    if (acc !== st.acc) return;
    st.editStamp++;
    st.sp.clear(); st.rc.clear();
    scheduleBase();
    renderStaleBanners();
    renderKoEdits();
    if (st.sweep) renderSweepTable();
  });

  // ------------------------------------------------------------- skeleton ----
  root.innerHTML = `
    <h1>Constraint-based analysis</h1>
    <p class="sub">Reaction knockouts, growth coupling and production envelopes for one GEM on one
    medium, solved with GLPK in this browser. Bound edits and knockouts made anywhere in this
    session apply to every solve. Methods and thresholds are on the
    <a href="methods.html#analysis">methods page</a>.</p>

    <div class="card">
      <div class="gem-toolbar">
        <div class="field">
          <label for="an-gem">GEM (1 of 35)</label>
          <select id="an-gem"><option value="">Choose a GEM…</option></select>
        </div>
        <div class="field">
          <label for="an-medium">Medium (${Object.keys(mediaLib).length} predefined)</label>
          <select id="an-medium"><option value="">Choose a medium…</option></select>
        </div>
        <div class="field picker an-picker">
          <label for="an-sub">Substrate (carbon source)</label>
          <input type="text" id="an-sub" role="combobox" aria-expanded="false" aria-controls="an-sub-listbox"
                 aria-autocomplete="list" autocomplete="off" placeholder="Optional: swaps in as carbon source" disabled>
          <ul class="listbox" id="an-sub-listbox" role="listbox" hidden></ul>
        </div>
        <div class="field picker an-picker">
          <label for="an-prod">Target product</label>
          <input type="text" id="an-prod" role="combobox" aria-expanded="false" aria-controls="an-prod-listbox"
                 aria-autocomplete="list" autocomplete="off" placeholder="Name or id in this GEM" disabled>
          <ul class="listbox" id="an-prod-listbox" role="listbox" hidden></ul>
        </div>
        <div class="field">
          <label for="an-frac">Growth fraction for the floor (%)</label>
          <input type="number" id="an-frac" min="50" max="99" step="1" value="90" style="width:90px">
        </div>
      </div>
      <div class="m2-swap" id="an-swapline"></div>
      <p class="status" id="an-medstatus" role="status" aria-live="polite"></p>
      <div id="an-koedits"></div>
    </div>

    <div class="an-verdict card" id="an-base" role="status" aria-live="polite">
      <p class="status">Choose a GEM and a medium above. The reference state (max growth, max
      product export, guaranteed product floor, coupling verdict) is solved as soon as both are set.</p>
    </div>

    <div class="card an-section" id="an-couple">
      <h2>Growth-coupling knockout search</h2>
      <p class="sub">Greedy search for a knockout set that makes the product's minimum export flux
      positive when biomass is held at or above the chosen fraction of the knockout strain's own
      maximum. Candidates are reactions carrying flux in a parsimonious zero-product optimum; the
      search applies the single knockout that most lowers decoupled growth and repeats. A failure
      to couple within K knockouts is reported as exactly that; the heuristic does not prove
      impossibility.</p>
      <div class="gem-toolbar">
        <div class="field">
          <label for="an-k">Max knockouts K</label>
          <input type="number" id="an-k" min="1" max="6" step="1" value="3" style="width:80px">
        </div>
        <div class="field">
          <label for="an-viab">Viability floor (% of reference max growth)</label>
          <input type="number" id="an-viab" min="1" max="90" step="1" value="10" style="width:90px">
        </div>
        <label class="mode" style="border-style:solid"><input type="checkbox" id="an-nogene">
          Allow knockouts of reactions without a gene rule</label>
        <button class="btn primary" id="an-couple-run" type="button" disabled>Search knockout set</button>
        <button class="btn small" id="an-couple-cancel" type="button" hidden>Cancel</button>
      </div>
      <div id="an-couple-prog" class="status" role="status" aria-live="polite"></div>
      <div id="an-couple-out"></div>
    </div>

    <div class="card an-section" id="an-envelope">
      <h2>Production envelope</h2>
      <p class="sub">Maximum and minimum product flux at each biomass level from 0 to max growth,
      under the current bounds. With a coupled knockout set from the search above, the knockout
      envelope overlays the reference one so the coupling is visible: a positive minimum at high
      biomass is the coupling.</p>
      <div class="gem-toolbar">
        <div class="field">
          <label for="an-envn">Biomass levels</label>
          <input type="number" id="an-envn" min="8" max="60" step="1" value="24" style="width:80px">
        </div>
        <button class="btn primary" id="an-env-run" type="button" disabled>Compute envelope</button>
        <button class="btn small" id="an-env-cancel" type="button" hidden>Cancel</button>
        <span class="status" id="an-env-prog" role="status" aria-live="polite"></span>
      </div>
      <div id="an-env-out"></div>
    </div>

    <div class="card an-section" id="an-sweep">
      <h2>Single-reaction knockout sweep</h2>
      <p class="sub">Every non-exchange reaction is closed in turn and the LP re-solved: max growth
      per knockout, and with a target product set, the knockout's max product export and its
      guaranteed floor at the growth fraction above. Each row is a fresh solve; the sweep runs
      incrementally and can be cancelled, keeping what it computed.</p>
      <div class="gem-toolbar">
        <button class="btn primary" id="an-sweep-run" type="button" disabled>Run knockout sweep</button>
        <button class="btn small" id="an-sweep-cancel" type="button" hidden>Cancel sweep</button>
        <label class="mode" style="border-style:solid"><input type="checkbox" id="an-sweep-filtered">
          Only reactions matching the table filter</label>
        <span class="status" id="an-sweep-scope"></span>
      </div>
      <div id="an-sweep-prog" class="status" role="status" aria-live="polite"></div>
      <div class="tablebar">
        <div class="field" style="flex:1 1 220px">
          <label for="an-sweep-filter">Filter (id, name, gene)</label>
          <input type="search" id="an-sweep-filter" placeholder="e.g. PGK or Glycolysis gene">
        </div>
        <div class="field">
          <label for="an-sweep-cls">Class</label>
          <select id="an-sweep-cls">
            <option value="">all classes</option>
            <option value="lethal">lethal</option>
            <option value="beneficial">raises product floor</option>
            <option value="costly">costly</option>
            <option value="neutral">neutral</option>
          </select>
        </div>
        <button class="btn small" id="an-sweep-csv" type="button" disabled>Export CSV</button>
        <button class="btn small" id="an-sweep-json" type="button" disabled>Export JSON</button>
      </div>
      <p class="count" id="an-sweep-count" role="status" aria-live="polite"></p>
      <div class="tablewrap"><table class="data" id="an-sweep-table">
        <thead><tr>
          <th scope="col"><button class="thsort" type="button" data-sort="id">Reaction</button></th>
          <th scope="col">Name</th>
          <th scope="col"><button class="thsort" type="button" data-sort="mu">KO max growth (h<sup>-1</sup>)</button></th>
          <th scope="col"><button class="thsort" type="button" data-sort="dmu">Change vs reference</button></th>
          <th scope="col"><button class="thsort" type="button" data-sort="productMax">KO max product (${UNIT})</button></th>
          <th scope="col"><button class="thsort" type="button" data-sort="floor">KO product floor (${UNIT})</button></th>
          <th scope="col"><button class="thsort" type="button" data-sort="cls">Class</button></th>
          <th scope="col"></th>
        </tr></thead>
        <tbody></tbody>
      </table></div>
      <div class="pager" style="margin-top:var(--s3)">
        <button class="btn small" id="an-sweep-prev" type="button">Previous</button>
        <span id="an-sweep-pageinfo"></span>
        <button class="btn small" id="an-sweep-next" type="button">Next</button>
      </div>

      <h3>Shadow prices and reduced costs (finite-difference estimates)</h3>
      <p class="sub">glpk.js does not expose LP duals, so these are estimated at the max-growth
      optimum by re-solving with the metabolite's mass-balance right-hand side moved by
      ${SHADOW_EPS}: one LP per metabolite. Reduced costs derive from those estimates over each
      reaction's own metabolites. At a degenerate optimum the estimate is a one-sided derivative;
      a metabolite whose perturbation is infeasible in both directions (a conserved pool) has no
      finite estimate and is reported so. Estimates within 10<sup>-6</sup> of zero should be read
      as zero.</p>
      <div class="an-duals">
        <div>
          <div class="field">
            <label for="an-sp-filter">Shadow prices: filter metabolites (id, name)</label>
            <input type="search" id="an-sp-filter" placeholder="e.g. atp or pyruvate">
          </div>
          <div class="tablebar">
            <button class="btn small" id="an-sp-run" type="button" disabled>Estimate</button>
            <button class="btn small" id="an-sp-csv" type="button" disabled>Export CSV</button>
            <span class="status" id="an-sp-prog" role="status" aria-live="polite"></span>
          </div>
          <div class="tablewrap" style="max-height:300px;overflow-y:auto"><table class="data" id="an-sp-table">
            <thead><tr><th scope="col">Metabolite</th><th scope="col">Name</th>
            <th scope="col">Shadow price est. (h<sup>-1</sup> per ${UNIT})</th></tr></thead>
            <tbody></tbody></table></div>
        </div>
        <div>
          <div class="field">
            <label for="an-rc-filter">Reduced costs: filter reactions (id, name)</label>
            <input type="search" id="an-rc-filter" placeholder="e.g. PYK or transport">
          </div>
          <div class="tablebar">
            <button class="btn small" id="an-rc-run" type="button" disabled>Estimate</button>
            <button class="btn small" id="an-rc-csv" type="button" disabled>Export CSV</button>
            <span class="status" id="an-rc-prog" role="status" aria-live="polite"></span>
          </div>
          <div class="tablewrap" style="max-height:300px;overflow-y:auto"><table class="data" id="an-rc-table">
            <thead><tr><th scope="col">Reaction</th><th scope="col">Name</th>
            <th scope="col">Reduced cost est. (h<sup>-1</sup> per ${UNIT})</th></tr></thead>
            <tbody></tbody></table></div>
        </div>
      </div>
    </div>`;

  const $ = (sel) => root.querySelector(sel);

  // ------------------------------------------------------------- controls ----
  const gemSel = $('#an-gem');
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
  const medSel = $('#an-medium');
  for (const label of Object.keys(mediaLib)) {
    const o = document.createElement('option');
    o.value = label; o.textContent = label;
    medSel.appendChild(o);
  }

  gemSel.addEventListener('change', async () => {
    const acc = gemSel.value;
    st.gem = null; st.acc = null; st.sub = null; st.prod = null;
    invalidateAll();
    $('#an-sub').value = ''; $('#an-prod').value = '';
    if (!acc) { refreshEnableState(); renderBase(); return; }
    try {
      $('#an-medstatus').textContent = `Loading GEM ${acc} (about 1 MB)…`;
      const gem = await loadGem(acc);
      st.gem = gem; st.acc = acc;
      ctx.setAccent(gem.species);
      $('#an-medstatus').textContent = `${gem.species} · ${acc} · ${fmt.format(gem.reactions.length)} reactions`;
      buildPickers();
      prefillEndpoints();
    } catch (e) {
      $('#an-medstatus').textContent = `Could not load GEM ${acc}: ${e.message}. Choose it again to retry.`;
    }
    refreshEnableState();
    renderSwapLine();
    renderKoEdits();
    scheduleBase();
  });

  medSel.addEventListener('change', () => {
    st.mediumLabel = medSel.value || null;
    invalidateAll();
    refreshEnableState();
    renderSwapLine();
    scheduleBase();
  });

  $('#an-frac').addEventListener('change', () => {
    const v = parseFloat($('#an-frac').value);
    if (!Number.isNaN(v) && v >= 50 && v <= 99) {
      st.frac = v / 100;
      invalidateAll();
      scheduleBase();
    }
  });

  // ------------------------------------------------------------- pickers ----
  // Lightweight combobox over the chosen GEM's metabolites (id + name).
  const metDisplayName = (m) => (m.name && m.name !== m.id ? m.name : '');

  function buildPickers() {
    const items = st.gem.metabolites.map(m => ({ mid: m.id, name: metDisplayName(m) }));
    makePicker('sub', items);
    makePicker('prod', items);
    $('#an-sub').disabled = false;
    $('#an-prod').disabled = false;
  }

  function makePicker(kind, items) {
    const input = $(`#an-${kind}`);
    const listbox = $(`#an-${kind}-listbox`);
    let options = [], active = -1;
    const close = () => { listbox.hidden = true; input.setAttribute('aria-expanded', 'false'); active = -1; };
    const open = () => { listbox.hidden = false; input.setAttribute('aria-expanded', 'true'); };
    const choose = (mid) => {
      st[kind === 'sub' ? 'sub' : 'prod'] = mid;
      const it = items.find(x => x.mid === mid);
      input.value = it && it.name ? `${it.name} (${mid})` : mid;
      close();
      invalidateAll();
      renderSwapLine();
      scheduleBase();
    };
    const render = (q) => {
      const query = q.trim().toLowerCase();
      if (!query) { close(); return; }
      const scored = [];
      for (const m of items) {
        const idL = m.mid.toLowerCase(), nmL = m.name.toLowerCase();
        let s = -1;
        if (idL === query || nmL === query) s = 0;
        else if (idL.startsWith(query) || nmL.startsWith(query)) s = 1;
        else if (idL.includes(query) || nmL.includes(query)) s = 2;
        if (s >= 0) scored.push([s, m]);
      }
      scored.sort((a, b) => a[0] - b[0] || a[1].mid.localeCompare(b[1].mid));
      const all = scored.map(x => x[1]);
      options = all.slice(0, 40);
      if (!all.length) {
        listbox.innerHTML = `<li class="mcap" role="presentation">No metabolite of ${esc(st.acc)} matches "${esc(q)}" (searching ${fmt.format(items.length)} ids and names).</li>`;
        open(); return;
      }
      const cap = all.length > options.length
        ? `<li class="mcap" role="presentation">Showing ${options.length} of ${fmt.format(all.length)} matches; keep typing to narrow.</li>` : '';
      listbox.innerHTML = options.map((m, i) => `
        <li id="an-${kind}-opt-${i}" role="option" aria-selected="false" data-mid="${esc(m.mid)}">
          <span class="mname-primary">${esc(m.name || m.mid)}</span>
          ${m.name ? `<span class="mid">${esc(m.mid)}</span>` : ''}
        </li>`).join('') + cap;
      listbox.querySelectorAll('li[role="option"]').forEach(li => {
        li.addEventListener('mousedown', (e) => { e.preventDefault(); choose(li.dataset.mid); });
      });
      open();
    };
    const setActive = (i) => {
      const lis = listbox.querySelectorAll('li[role="option"]');
      if (!lis.length) return;
      active = (i + lis.length) % lis.length;
      lis.forEach((li, j) => li.setAttribute('aria-selected', String(j === active)));
      input.setAttribute('aria-activedescendant', `an-${kind}-opt-${active}`);
      lis[active].scrollIntoView({ block: 'nearest' });
    };
    input.oninput = () => { st[kind === 'sub' ? 'sub' : 'prod'] = null; render(input.value); };
    input.onblur = () => setTimeout(close, 120);
    input.onkeydown = (e) => {
      if (listbox.hidden && (e.key === 'ArrowDown' || e.key === 'ArrowUp')) { render(input.value); return; }
      if (e.key === 'ArrowDown') { e.preventDefault(); setActive(active + 1); }
      else if (e.key === 'ArrowUp') { e.preventDefault(); setActive(active - 1); }
      else if (e.key === 'Enter') {
        if (!listbox.hidden && active >= 0 && options[active]) { e.preventDefault(); choose(options[active].mid); }
      } else if (e.key === 'Escape') close();
    };
  }

  function prefillEndpoints() {
    const ep = ctx.getEndpoints();
    const has = (mid) => mid && st.gem.metabolites.some(m => m.id === mid);
    if (has(ep.sub)) {
      st.sub = ep.sub;
      const m = st.gem.metabolites.find(x => x.id === ep.sub);
      const nm = metDisplayName(m);
      $('#an-sub').value = nm ? `${nm} (${ep.sub})` : ep.sub;
    }
    if (has(ep.prod)) {
      st.prod = ep.prod;
      const m = st.gem.metabolites.find(x => x.id === ep.prod);
      const nm = metDisplayName(m);
      $('#an-prod').value = nm ? `${nm} (${ep.prod})` : ep.prod;
    }
  }

  // ------------------------------------------------------ medium assembly ----
  function workingBounds() {
    const def = mediaLib[st.mediumLabel];
    if (!def) return null;
    const out = {};
    for (const [ex, lb] of Object.entries(def.components || {})) if (!ex.startsWith('_')) out[ex] = lb;
    for (const [ex, lb] of Object.entries(def.supplements || {})) if (!ex.startsWith('_') && !(ex in out)) out[ex] = lb;
    const info = { carbonEx: def.carbon_exchange || null, cap: def.carbon_cap ?? null, swapped: false, subEx: null, closedEx: null };
    if (st.swapCarbon && st.sub && st.gem) {
      const base = st.sub.replace(/_[a-z]+$/, '');
      const exId = `EX_${base}_e`;
      if (st.gem.reactions.some(r => r.ex && r.id === exId)) {
        if (info.carbonEx && info.carbonEx !== exId) { out[info.carbonEx] = 0; info.closedEx = info.carbonEx; }
        out[exId] = info.cap ?? -10;
        info.swapped = true; info.subEx = exId;
      }
    }
    return { bounds: out, info, def };
  }

  function substrateExchangeForYield(info) {
    return (info && (info.subEx || info.carbonEx)) || null;
  }

  function renderSwapLine() {
    const host = $('#an-swapline');
    if (!st.mediumLabel || !st.gem) { host.innerHTML = ''; return; }
    const w = workingBounds();
    let line;
    if (!st.sub) {
      line = `Medium as defined: carbon source ${w.info.carbonEx ? `<span class="mono">${esc(w.info.carbonEx)}</span> at lb ${w.info.cap ?? 'not recorded'}` : 'not recorded'}. Pick a substrate to swap it in.`;
    } else if (!st.swapCarbon) {
      line = `The substrate is not swapped in; the medium's carbon source ${w.info.carbonEx ? `<span class="mono">${esc(w.info.carbonEx)}</span>` : ''} stays open.`;
    } else if (!w.info.swapped) {
      line = `${esc(st.acc)} has no exchange for <span class="mono">${esc(st.sub)}</span>; the medium runs as defined.`;
    } else if (w.info.subEx === w.info.carbonEx) {
      line = `The substrate already is the medium's carbon source (<span class="mono">${esc(w.info.subEx)}</span> at lb ${w.info.cap ?? 'not recorded'}).`;
    } else {
      line = `${w.info.closedEx ? `<span class="mono">${esc(w.info.closedEx)}</span> closes and ` : ''}<span class="mono">${esc(w.info.subEx)}</span> opens at lb ${w.info.cap ?? -10}.`;
    }
    host.innerHTML = `
      <label class="mode" style="border-style:solid"><input type="checkbox" id="an-swap-cb" ${st.swapCarbon ? 'checked' : ''}>
        Substrate replaces the medium's carbon source</label>
      <span class="status">${line}</span>`;
    host.querySelector('#an-swap-cb').addEventListener('change', (e) => {
      st.swapCarbon = e.target.checked;
      invalidateAll();
      renderSwapLine();
      scheduleBase();
    });
  }

  // Knockouts held in the session's bound edits (lb = ub = 0).
  function sessionKOs() {
    return st.acc ? listEdits(st.acc).filter(e => e.lb === 0 && e.ub === 0).map(e => e.rid) : [];
  }
  function renderKoEdits() {
    const host = $('#an-koedits');
    if (!st.acc) { host.innerHTML = ''; return; }
    const edits = listEdits(st.acc);
    const kos = edits.filter(e => e.lb === 0 && e.ub === 0);
    if (!edits.length) { host.innerHTML = ''; return; }
    host.innerHTML = `<p class="count">${edits.length} bound edit${edits.length > 1 ? 's' : ''} in this session
      (${kos.length} knockout${kos.length === 1 ? '' : 's'}), applied to every solve:</p>
      <div class="chiprow">${kos.map(k => `
        <span class="an-kochip"><span class="mono">${esc(k.rid)}</span>
          <button class="btn small" type="button" data-restore="${esc(k.rid)}" aria-label="Restore ${esc(k.rid)}">Restore</button>
        </span>`).join('')}</div>`;
    host.querySelectorAll('button[data-restore]').forEach(b =>
      b.addEventListener('click', () => clearEdit(st.acc, b.dataset.restore)));
  }

  // -------------------------------------------------------------- sessions ----
  function ready() { return !!(st.gem && st.mediumLabel); }
  async function newSession(extraKos = []) {
    const w = workingBounds();
    const target = st.prod ? analysisTarget(st.gem, st.prod) : null;
    const S = await makeSession(st.gem, st.acc, w.bounds, target, { kos: extraKos });
    return { S, w, target };
  }

  function refreshEnableState() {
    const ok = ready();
    const okP = ok && !!st.prod;
    $('#an-sweep-run').disabled = !ok || !!st.busy;
    $('#an-couple-run').disabled = !okP || !!st.busy;
    $('#an-env-run').disabled = !okP || !!st.busy;
    $('#an-sp-run').disabled = !ok || !!st.busy;
    $('#an-rc-run').disabled = !ok || !!st.busy;
  }

  function invalidateAll() {
    st.base = null;
    st.sp.clear(); st.rc.clear();
    renderStaleBanners();
  }

  function renderStaleBanners() {
    for (const [key, obj] of [['#an-couple-out', st.search], ['#an-env-out', st.env], ['#an-sweep-prog', st.sweep]]) {
      const el = $(key);
      if (!el || !obj) continue;
      let b = el.querySelector('.an-stale');
      const stale = obj.stamp !== staleStamp();
      if (stale && !b) {
        b = document.createElement('p');
        b.className = 'termination capped an-stale';
        b.textContent = 'The GEM, medium, target or bounds changed after this result was computed; re-run to refresh it.';
        el.prepend(b);
      } else if (!stale && b) b.remove();
    }
  }
  function staleStamp() {
    return [st.acc, st.mediumLabel, st.sub, st.prod, st.swapCarbon, st.frac, st.editStamp].join('|');
  }

  // ------------------------------------------------------------ base state ----
  let baseTimer = null;
  function scheduleBase() {
    refreshEnableState();
    if (baseTimer) clearTimeout(baseTimer);
    baseTimer = setTimeout(runBase, 250);
  }

  async function runBase() {
    const token = ++st.baseToken;
    const el = $('#an-base');
    if (!ready()) { renderBase(); return; }
    el.innerHTML = `<span class="status">Solving the reference state of ${esc(st.acc)} on ${esc(st.mediumLabel)} (${st.prod ? 'up to 3' : '1'} LPs)…</span>`;
    try {
      const { S, w, target } = await newSession();
      if (st.prod && !target) {
        if (token !== st.baseToken) return;
        st.base = { noTarget: true };
        renderBase();
        return;
      }
      const b = await baseState(S, st.frac);
      if (token !== st.baseToken) return;
      st.base = { ...b, target, w, stamp: staleStamp() };
      renderBase();
      renderStaleBanners();
    } catch (e) {
      if (token !== st.baseToken) return;
      el.innerHTML = `<p class="status error">Reference solve failed: ${esc(e.message)}. Change a setting to retry.</p>`;
    }
    refreshEnableState();
  }

  function targetLine(target) {
    if (!target) return '';
    if (target.kind === 'exchange') {
      return `Product measured as export flux <span class="mono">${esc(target.id)}</span>${target.mid && !target.mid.endsWith('_e') ? ` (the exchange of the extracellular form of <span class="mono">${esc(target.mid)}</span>)` : ''}.`;
    }
    return `No exchange for this product in ${esc(st.acc)}; it is measured on an added demand <span class="mono">${esc(target.id)}</span>, the only sink, so its minimum is the guaranteed net production.`;
  }

  function renderBase() {
    const el = $('#an-base');
    if (!ready()) {
      el.innerHTML = `<p class="status">Choose a GEM and a medium above. The reference state (max growth, max
      product export, guaranteed product floor, coupling verdict) is solved as soon as both are set.</p>`;
      return;
    }
    const b = st.base;
    if (!b) return;   // solve in flight; runBase already rendered the loading line
    if (b.noTarget) {
      el.innerHTML = `<p class="status error">${esc(st.acc)} does not contain the chosen product metabolite; pick one from the product list (it searches this GEM only).</p>`;
      return;
    }
    const kos = sessionKOs();
    const parts = [];
    parts.push(chipHTML(b.mu != null && b.mu > ZERO_MU ? 'ok' : 'bad',
      b.mu == null ? `max growth: ${esc(statusName(b.muStatus))}` : `max growth ${fnum(b.mu, 3)} h<sup>-1</sup>`));
    if (b.target !== undefined && st.prod) {
      parts.push(chipHTML(b.productMax != null && b.productMax > FLUX_TOL ? 'ok' : 'bad',
        b.productMax == null ? `max product: ${esc(statusName(b.productStatus))}` : `max product ${fnum(b.productMax, 3)} ${UNIT}`));
      const fl = b.floor;
      parts.push(chipHTML(b.coupled ? 'ok' : (b.coupled === false ? 'bad' : 'na'),
        fl == null
          ? `product floor at ≥${Math.round(st.frac * 100)}% growth: not computed${b.mu != null && b.mu <= ZERO_MU ? ' (no growth)' : ''}`
          : `growth-coupled: ${b.coupled ? 'yes' : 'no'} (floor ${fnum(fl, 4)} ${UNIT} at ≥${Math.round(st.frac * 100)}% of max growth)`));
    }
    el.innerHTML = `
      <div class="an-chips">${parts.join('')}</div>
      <p class="status">Reference state of ${esc(st.acc)} on ${esc(st.mediumLabel)}${kos.length ? ` with ${kos.length} session knockout${kos.length === 1 ? '' : 's'}` : ''}${st.base.w && st.base.w.info.swapped ? ` · substrate swapped in as carbon source (<span class="mono">${esc(st.base.w.info.subEx)}</span>)` : ''} · ${b.solves} LP${b.solves === 1 ? '' : 's'} solved in this browser.
      ${st.prod ? targetLine(b.target) : 'Pick a target product to add product-side results (max export, floor, coupling).'}</p>
      ${st.prod && b.target && b.target.kind === 'exchange' && st.base.w && (st.base.w.bounds[b.target.id] ?? 0) < 0
        ? `<p class="termination capped">The medium allows uptake of the product (<span class="mono">${esc(b.target.id)}</span> lb ${st.base.w.bounds[b.target.id]}); a negative floor means net uptake, not production.</p>` : ''}`;
  }

  function chipHTML(cls, html) { return `<span class="feas-chip ${cls}">${html}</span>`; }

  // ---------------------------------------------------------- long-run gate ----
  function startRun(name) {
    st.busy = name; st.cancelFlag = false;
    refreshEnableState();
  }
  function endRun() {
    st.busy = null; st.cancelFlag = false;
    refreshEnableState();
  }

  // -------------------------------------------------------------- KO sweep ----
  const sweepUI = {
    page: 0, pageSize: 50, sortKey: 'id', sortDir: 1, filter: '', cls: '',
  };

  $('#an-sweep-filter').addEventListener('input', (e) => { sweepUI.filter = e.target.value.trim().toLowerCase(); sweepUI.page = 0; renderSweepTable(); });
  $('#an-sweep-cls').addEventListener('change', (e) => { sweepUI.cls = e.target.value; sweepUI.page = 0; renderSweepTable(); });
  $('#an-sweep-prev').addEventListener('click', () => { if (sweepUI.page > 0) { sweepUI.page--; renderSweepTable(); } });
  $('#an-sweep-next').addEventListener('click', () => { sweepUI.page++; renderSweepTable(); });
  root.querySelectorAll('#an-sweep-table .thsort').forEach(btn => btn.addEventListener('click', () => {
    const k = btn.dataset.sort;
    if (sweepUI.sortKey === k) sweepUI.sortDir *= -1; else { sweepUI.sortKey = k; sweepUI.sortDir = 1; }
    renderSweepTable();
  }));
  $('#an-sweep-cancel').addEventListener('click', () => { st.cancelFlag = true; });
  $('#an-sweep-csv').addEventListener('click', () => exportSweep('csv'));
  $('#an-sweep-json').addEventListener('click', () => exportSweep('json'));

  function sweepScopeNote() {
    if (!st.gem) return '';
    const { rxns, excluded, total } = sweepScope(st.gem, st.gem.stats.biomass_id);
    return `Scope: ${fmt.format(rxns.length)} of ${fmt.format(total)} reactions (${fmt.format(excluded)} exchanges and the biomass reaction excluded).`;
  }

  $('#an-sweep-run').addEventListener('click', async () => {
    if (!ready() || st.busy) return;
    startRun('sweep');
    const prog = $('#an-sweep-prog');
    const btnCancel = $('#an-sweep-cancel');
    btnCancel.hidden = false;
    const t0 = performance.now();
    try {
      const { S, target } = await newSession();
      const refs = await baseState(S, st.frac);
      refs.frac = st.frac;
      const onlyFiltered = $('#an-sweep-filtered').checked && sweepUI.filter;
      const scope = sweepScope(st.gem, st.gem.stats.biomass_id);
      let restrictTo = null;
      if (onlyFiltered) {
        restrictTo = new Set(scope.rxns.filter(r => rowMatchesFilter(r, sweepUI.filter)).map(r => r.id));
        if (!restrictTo.size) {
          prog.textContent = `No reaction in the sweep scope matches the filter "${sweepUI.filter}"; nothing to sweep.`;
          endRun(); btnCancel.hidden = true; return;
        }
      }
      st.sweep = {
        rows: [], refs, stamp: staleStamp(), cancelled: false,
        scope: restrictTo ? restrictTo.size : scope.rxns.length,
        fullScope: scope.rxns.length, excluded: scope.excluded, total: scope.total,
        restricted: !!restrictTo, filterUsed: onlyFiltered ? sweepUI.filter : null,
        target: !!target,
      };
      renderSweepTable();
      let lastRender = 0;
      const skip = restrictTo ? (rid) => !restrictTo.has(rid) : null;
      const res = await koSweepFiltered(S, refs, skip, {
        shouldStop: () => st.cancelFlag,
        cbRow: (row) => {
          st.sweep.rows.push(row);
          const now = performance.now();
          const rate = (now - t0) / Math.max(S.solves, 1);
          prog.innerHTML = `Sweeping: ${st.sweep.rows.length} of ${st.sweep.scope} knockouts solved · ${S.solves} LPs ·
            ${rate.toFixed(0)} ms per LP · about ${etaText((st.sweep.scope - st.sweep.rows.length) / Math.max(st.sweep.rows.length, 1) * (now - t0))} remaining.`;
          if (now - lastRender > 600) { lastRender = now; renderSweepTable(); }
        },
      });
      st.sweep.cancelled = res.cancelled;
      st.sweep.solves = res.solves;
      const secs = ((performance.now() - t0) / 1000).toFixed(0);
      prog.innerHTML = res.cancelled
        ? `Sweep cancelled: ${st.sweep.rows.length} of ${st.sweep.scope} knockouts solved (${res.solves} LPs, ${secs} s). Computed rows stay below and export as a partial sweep.`
        : `Sweep complete: ${st.sweep.rows.length} of ${st.sweep.scope} knockouts solved (${res.solves} LPs, ${secs} s).`;
      renderSweepTable();
    } catch (e) {
      prog.textContent = `Sweep failed: ${e.message}. Computed rows, if any, stay below.`;
    } finally {
      btnCancel.hidden = true;
      endRun();
      $('#an-sweep-csv').disabled = !st.sweep || !st.sweep.rows.length;
      $('#an-sweep-json').disabled = !st.sweep || !st.sweep.rows.length;
    }
  });

  // Wrap the engine sweep so a filter restriction can skip rows without
  // spending solves on them.
  async function koSweepFiltered(S, refs, skip, opts) {
    if (!skip) return koSweep(S, refs, opts);
    // Restricted sweep: run the engine on a session view whose GEM only lists
    // the reactions in scope, so skipped rows cost no solves. The LP itself is
    // shared and untouched, so results are identical to a full-model sweep of
    // the same reactions.
    const reduced = { ...S.gem, reactions: S.gem.reactions.filter(r => r.ex || r.id === S.biomass || !skip(r.id)) };
    const S2 = Object.create(S);
    S2.gem = reduced;
    return koSweep(S2, refs, opts);
  }

  function etaText(ms) {
    if (!Number.isFinite(ms) || ms < 0) return 'unknown time';
    const s = Math.round(ms / 1000);
    if (s < 90) return `${s} s`;
    return `${Math.round(s / 60)} min`;
  }

  function rowMatchesFilter(r, q) {
    return r.id.toLowerCase().includes(q) || (r.name || '').toLowerCase().includes(q)
      || (r.genes || []).some(g => g.toLowerCase().includes(q));
  }

  function filteredSweepRows() {
    if (!st.sweep) return [];
    let rows = st.sweep.rows;
    if (sweepUI.filter) rows = rows.filter(r => rowMatchesFilter(r, sweepUI.filter));
    if (sweepUI.cls) rows = rows.filter(r => r.cls === sweepUI.cls);
    const k = sweepUI.sortKey, dir = sweepUI.sortDir;
    const val = (r) => {
      if (k === 'dmu') return (r.mu == null || st.sweep.refs.mu == null) ? -Infinity : (r.mu - st.sweep.refs.mu);
      if (k === 'id' || k === 'cls') return r[k] || '';
      return r[k] == null ? -Infinity : r[k];
    };
    return [...rows].sort((a, b) => {
      const va = val(a), vb = val(b);
      const c = (typeof va === 'string') ? va.localeCompare(vb) : (va - vb);
      return dir * (c || a.id.localeCompare(b.id));
    });
  }

  function renderSweepTable() {
    $('#an-sweep-scope').textContent = sweepScopeNote();
    const tbody = $('#an-sweep-table tbody');
    if (!st.sweep) {
      tbody.innerHTML = '';
      $('#an-sweep-count').textContent = st.gem
        ? 'No sweep run yet for this GEM and medium.'
        : 'Choose a GEM and a medium, then run the sweep.';
      $('#an-sweep-pageinfo').textContent = '';
      $('#an-sweep-prev').disabled = $('#an-sweep-next').disabled = true;
      return;
    }
    const list = filteredSweepRows();
    const pages = Math.max(1, Math.ceil(list.length / sweepUI.pageSize));
    sweepUI.page = Math.min(sweepUI.page, pages - 1);
    const slice = list.slice(sweepUI.page * sweepUI.pageSize, (sweepUI.page + 1) * sweepUI.pageSize);
    const refs = st.sweep.refs;
    const kos = new Set(sessionKOs());
    $('#an-sweep-count').textContent =
      `Showing ${fmt.format(slice.length)} of ${fmt.format(list.length)} ` +
      (sweepUI.filter || sweepUI.cls ? `matching knockouts (filtered from ${fmt.format(st.sweep.rows.length)} computed)` : 'computed knockouts') +
      ` · ${fmt.format(st.sweep.rows.length)} of ${fmt.format(st.sweep.scope)} in scope computed` +
      (st.sweep.restricted ? ` (sweep restricted to filter "${st.sweep.filterUsed}"; full scope ${fmt.format(st.sweep.fullScope)})` : '') +
      ` · reference max growth ${fnum(refs.mu, 4)} h^-1` +
      (st.sweep.target ? `, reference floor ${fnum(refs.floor, 4)}` : ', no product target set');
    $('#an-sweep-pageinfo').textContent = `Page ${sweepUI.page + 1} of ${pages}`;
    $('#an-sweep-prev').disabled = sweepUI.page === 0;
    $('#an-sweep-next').disabled = sweepUI.page >= pages - 1;

    tbody.innerHTML = slice.map(r => {
      const dmu = (r.mu == null || refs.mu == null || refs.mu <= ZERO_MU) ? null : (r.mu - refs.mu) / refs.mu * 100;
      const isKO = kos.has(r.id);
      return `<tr>
        <td class="mono">${esc(r.id)}${r.alreadyClosed ? ' <span class="an-note">bounds already 0</span>' : ''}</td>
        <td>${esc(r.name)}</td>
        <td class="mono">${r.lethalNoSteadyState ? 'no steady state' : fnum(r.mu)}</td>
        <td class="mono">${dmu == null ? '' : (dmu >= 0 ? '+' : '') + dmu.toFixed(1) + '%'}</td>
        <td class="mono">${st.sweep.target ? fnum(r.productMax, 3) : ''}</td>
        <td class="mono">${st.sweep.target ? fnum(r.floor, 4) : ''}</td>
        <td><span class="an-cls ${esc(r.cls)}">${esc(CLS_LABEL[r.cls] || r.cls)}</span></td>
        <td><button class="btn small" type="button" data-ko="${esc(r.id)}">${isKO ? 'Restore' : 'KO in session'}</button></td>
      </tr>`;
    }).join('');
    tbody.querySelectorAll('button[data-ko]').forEach(b => b.addEventListener('click', () => {
      const rid = b.dataset.ko;
      if (new Set(sessionKOs()).has(rid)) clearEdit(st.acc, rid);
      else setEdit(st.acc, rid, 0, 0);
    }));
  }

  function exportSweep(kind) {
    if (!st.sweep || !st.sweep.rows.length) return;
    const refs = st.sweep.refs;
    const meta = {
      gem: st.acc, species: st.gem.species, medium: st.mediumLabel,
      substrate: st.sub, product: st.prod, growth_fraction: st.frac,
      reference_mu: refs.mu, reference_product_max: refs.productMax ?? null, reference_floor: refs.floor ?? null,
      knockouts_computed: st.sweep.rows.length, knockouts_in_scope: st.sweep.scope,
      scope_note: `${st.sweep.excluded} exchanges and the biomass reaction excluded from ${st.sweep.total} model reactions` +
        (st.sweep.restricted ? `; sweep restricted to filter "${st.sweep.filterUsed}" (full scope ${st.sweep.fullScope})` : ''),
      cancelled: st.sweep.cancelled,
      classification: `lethal: mu <= ${ZERO_MU}; costly: mu < ${COSTLY_FRAC} x reference; raises product floor: floor > reference floor + ${FLUX_TOL}`,
      units: 'mu in 1/h; product fluxes in mmol/gDW/h',
      session_bound_edits: listEdits(st.acc),
    };
    const rows = st.sweep.rows.map(r => ({
      reaction: r.id, name: r.name, genes: (r.genes || []).join(';'),
      ko_mu: r.lethalNoSteadyState ? 'no_steady_state' : r.mu,
      ko_product_max: r.productMax, ko_product_floor: r.floor,
      class: r.cls, bounds_already_zero: !!r.alreadyClosed,
    }));
    const base = `ko_sweep_${st.acc}_${(st.prod || 'no_product')}`;
    if (kind === 'json') {
      downloadBlob(JSON.stringify({ meta, rows }, null, 1), `${base}.json`, 'application/json');
    } else {
      const head = 'reaction,name,genes,ko_mu,ko_product_max,ko_product_floor,class,bounds_already_zero';
      const csv = [
        `# knockout sweep: ${meta.knockouts_computed} of ${meta.knockouts_in_scope} in scope${meta.cancelled ? ' (cancelled early)' : ''}; ${meta.scope_note}`,
        `# GEM ${meta.gem} on ${meta.medium}; reference mu ${meta.reference_mu}; product ${meta.product ?? 'none'}; ${meta.classification}`,
        head,
        ...rows.map(r => [r.reaction, r.name, r.genes, r.ko_mu, r.ko_product_max, r.ko_product_floor, r.class, r.bounds_already_zero].map(csvEscape).join(',')),
      ].join('\n');
      downloadBlob(csv, `${base}.csv`, 'text/csv');
    }
  }

  // ------------------------------------------------- shadow prices / costs ----
  $('#an-sp-run').addEventListener('click', () => runDuals('sp'));
  $('#an-rc-run').addEventListener('click', () => runDuals('rc'));
  $('#an-sp-filter').addEventListener('input', renderDualTables);
  $('#an-rc-filter').addEventListener('input', renderDualTables);
  $('#an-sp-csv').addEventListener('click', () => exportDuals('sp'));
  $('#an-rc-csv').addEventListener('click', () => exportDuals('rc'));

  function dualCandidates(kind) {
    const q = $(kind === 'sp' ? '#an-sp-filter' : '#an-rc-filter').value.trim().toLowerCase();
    if (kind === 'sp') {
      const all = st.gem ? st.gem.metabolites : [];
      return q ? all.filter(m => m.id.toLowerCase().includes(q) || (m.name || '').toLowerCase().includes(q)) : all;
    }
    const all = st.gem ? st.gem.reactions.filter(r => !r.ex) : [];
    return q ? all.filter(r => r.id.toLowerCase().includes(q) || (r.name || '').toLowerCase().includes(q)) : all;
  }

  async function runDuals(kind) {
    if (!ready() || st.busy) return;
    const cap = kind === 'sp' ? SP_BATCH : RC_BATCH;
    const cand = dualCandidates(kind).filter(x => !(kind === 'sp' ? st.sp : st.rc).has(x.id));
    const batch = cand.slice(0, cap).map(x => x.id);
    const prog = $(kind === 'sp' ? '#an-sp-prog' : '#an-rc-prog');
    if (!batch.length) {
      prog.textContent = cand.length ? '' : 'Every matching entry is already estimated (or nothing matches).';
      return;
    }
    startRun(kind);
    try {
      const { S } = await newSession();
      const opts = {
        shouldStop: () => st.cancelFlag,
        onProgress: (d, n) => { prog.textContent = `Estimating ${d} of ${n} (one LP each)…`; },
      };
      const matching = dualCandidates(kind).length;
      const left = cand.length - batch.length;
      if (kind === 'sp') {
        const res = await shadowPrices(S, batch, opts);
        if (!res.ok) { prog.textContent = `Reference LP not optimal (${res.statusText}); no estimates.`; return; }
        for (const [mid, rec] of res.values) st.sp.set(mid, rec);
      } else {
        const res = await reducedCosts(S, batch, opts);
        if (!res.ok) { prog.textContent = `Reference LP not optimal (${res.statusText}); no estimates.`; return; }
        for (const [rid, rec] of res.values) st.rc.set(rid, rec);
      }
      prog.textContent = `Estimated ${batch.length} this run; ${matching - left} of ${matching} matching ${kind === 'sp' ? 'metabolites' : 'reactions'} now estimated` +
        (left > 0 ? ` (${left} remain; run again for the next ${Math.min(cap, left)})` : '') + '.';
      renderDualTables();
    } catch (e) {
      prog.textContent = `Estimation failed: ${e.message}.`;
    } finally {
      endRun();
    }
  }

  function renderDualTables() {
    if (!st.gem) return;
    const spBody = $('#an-sp-table tbody');
    const rcBody = $('#an-rc-table tbody');
    const spCand = dualCandidates('sp');
    const rcCand = dualCandidates('rc');
    const spAll = spCand.filter(m => st.sp.has(m.id));
    const rcAll = rcCand.filter(r => st.rc.has(r.id));
    const spRows = spAll.slice(0, 400);
    const rcRows = rcAll.slice(0, 400);
    const capNote = (shown, all) => shown < all
      ? `<tr><td colspan="3" class="an-note">Showing ${shown} of ${all} estimated matches; narrow the filter to see the rest (export includes all).</td></tr>` : '';
    spBody.innerHTML = spRows.map(m => {
      const rec = st.sp.get(m.id);
      return `<tr><td class="mono">${esc(m.id)}</td><td>${esc(metDisplayName(m))}</td>
        <td class="mono">${rec.y == null ? `<span class="an-note">${esc(rec.note || 'not computed')}</span>` : rec.y.toExponential(3) + (rec.side === 'left' ? ' <span class="an-note">left derivative</span>' : '')}</td></tr>`;
    }).join('') + capNote(spRows.length, spAll.length)
      || `<tr><td colspan="3" class="an-note">${spCand.length ? `No estimates yet for the ${fmt.format(spCand.length)} matching metabolites; Estimate runs up to ${SP_BATCH} at a time.` : 'No metabolite matches the filter.'}</td></tr>`;
    rcBody.innerHTML = rcRows.map(r => {
      const rec = st.rc.get(r.id);
      return `<tr><td class="mono">${esc(r.id)}</td><td>${esc(r.name || '')}</td>
        <td class="mono">${rec.rc == null ? `<span class="an-note">${esc(rec.note || 'not computed')}</span>` : rec.rc.toExponential(3)}</td></tr>`;
    }).join('') + capNote(rcRows.length, rcAll.length)
      || `<tr><td colspan="3" class="an-note">${rcCand.length ? `No estimates yet for the ${fmt.format(rcCand.length)} matching reactions; Estimate runs up to ${RC_BATCH} at a time.` : 'No reaction matches the filter.'}</td></tr>`;
    $('#an-sp-csv').disabled = !st.sp.size;
    $('#an-rc-csv').disabled = !st.rc.size;
  }

  function exportDuals(kind) {
    const isSp = kind === 'sp';
    const store = isSp ? st.sp : st.rc;
    if (!store.size) return;
    const denom = isSp ? st.gem.metabolites.length : st.gem.reactions.filter(r => !r.ex).length;
    const head = isSp ? 'metabolite,shadow_price_estimate,side_or_note' : 'reaction,reduced_cost_estimate,note';
    const lines = [
      `# ${isSp ? 'shadow price' : 'reduced cost'} finite-difference estimates at the max-growth optimum; eps ${SHADOW_EPS}; GEM ${st.acc} on ${st.mediumLabel}; ${store.size} of ${denom} ${isSp ? 'metabolites' : 'non-exchange reactions'} estimated`,
      head,
      ...[...store.entries()].map(([id, rec]) => [
        id, isSp ? (rec.y ?? '') : (rec.rc ?? ''), isSp ? (rec.note || rec.side || '') : (rec.note || ''),
      ].map(csvEscape).join(',')),
    ];
    downloadBlob(lines.join('\n'), `${isSp ? 'shadow_prices' : 'reduced_costs'}_${st.acc}.csv`, 'text/csv');
  }

  // -------------------------------------------------------- coupling search ----
  $('#an-couple-cancel').addEventListener('click', () => { st.cancelFlag = true; });
  $('#an-couple-run').addEventListener('click', async () => {
    if (!ready() || !st.prod || st.busy) return;
    startRun('couple');
    const prog = $('#an-couple-prog');
    const out = $('#an-couple-out');
    const btnCancel = $('#an-couple-cancel');
    btnCancel.hidden = false;
    out.innerHTML = '';
    const t0 = performance.now();
    const K = clampInt($('#an-k').value, 1, 6, 3);
    const viab = clampInt($('#an-viab').value, 1, 90, 10) / 100;
    const iterLines = [];
    try {
      const { S, w, target } = await newSession();
      if (!target) {
        out.innerHTML = `<p class="status error">${esc(st.acc)} does not contain the chosen product; pick one from the product list.</p>`;
        return;
      }
      const muRefSolve = await baseState(S, st.frac);
      const res = await couplingSearch(S, {
        frac: st.frac, K, requireGene: !$('#an-nogene').checked,
        minMu: muRefSolve.mu != null ? viab * muRefSolve.mu : undefined,
        shouldStop: () => st.cancelFlag,
        onCandidate: (t, n, it) => { prog.textContent = `Iteration ${it}: scoring candidate knockout ${t} of ${n} (2 LPs each)…`; },
        onIteration: (it) => {
          iterLines.push(`Iteration ${it.iteration}: max growth ${fnum(it.mu, 4)}, decoupled growth ${fnum(it.mu0, 4)} (ratio ${fnum(it.ratio, 3)}), floor ${fnum(it.floor, 4)}; ${it.tested} of ${it.candidates} candidates viable, ${it.skippedLethal} below the viability floor` + (it.best ? `; best: ${it.best.id} (ratio ${fnum(it.best.ratio, 3)})` : '') + '.');
        },
      });
      const secs = ((performance.now() - t0) / 1000).toFixed(0);
      prog.textContent = '';
      st.search = { res, stamp: staleStamp(), w, target, secs };
      renderCoupleResult(S, res, w, target, secs, iterLines);
    } catch (e) {
      out.innerHTML = `<p class="status error">Coupling search failed: ${esc(e.message)}.</p>`;
    } finally {
      btnCancel.hidden = true;
      endRun();
    }
  });

  function clampInt(v, lo, hi, dflt) {
    const n = parseInt(v, 10);
    return Number.isNaN(n) ? dflt : Math.min(hi, Math.max(lo, n));
  }

  function renderCoupleResult(S, res, w, target, secs, iterLines) {
    const out = $('#an-couple-out');
    const coupled = res.verdict === 'coupled';
    const verdictChip = chipHTML(coupled ? 'ok' : 'bad',
      res.verdict === 'coupled' ? `growth-coupled: yes`
        : res.verdict === 'cancelled' ? 'search cancelled'
          : res.verdict === 'no-growth' ? 'no growth on this medium'
            : res.verdict === 'no-production' ? 'product not producible'
              : 'growth-coupled: no (within this search)');
    let koTable = '';
    if (res.koSet.length) {
      koTable = `
        <p class="count">${res.koSet.length} knockout${res.koSet.length === 1 ? '' : 's'} applied by the search (of at most ${res.K}):</p>
        <div class="tablewrap" style="max-width:760px"><table class="data">
          <thead><tr><th scope="col">Reaction</th><th scope="col">Name</th><th scope="col">Genes (GPR)</th></tr></thead>
          <tbody>${res.koSet.map(k => `<tr>
            <td class="mono">${esc(k.id)}</td><td>${esc(k.name)}</td>
            <td class="mono" style="overflow-wrap:anywhere">${esc(k.gpr || (k.genes || []).join(' ') || 'no gene rule')}</td>
          </tr>`).join('')}</tbody>
        </table></div>`;
    }
    let numbers = '';
    if (coupled) {
      const subEx = substrateExchangeForYield(w.info);
      const y = floorYield(S, res.floorVars, res.floor, subEx, subEx ? subEx.replace(/^EX_/, '') : null, st.prod);
      numbers = `
        <div class="an-chips">
          ${chipHTML('ok', `max growth with knockouts ${fnum(res.mu, 4)} h<sup>-1</sup>`)}
          ${chipHTML('ok', `guaranteed product ≥ ${fnum(res.floor, 4)} ${UNIT} at ≥${Math.round(res.frac * 100)}% of max growth`)}
          ${chipHTML('na', `max product ${fnum(res.productMax, 3)} ${UNIT}`)}
        </div>
        <p class="status">Yield at the floor solution: ${y.mmol == null ? esc(y.note)
          : `${y.mmol.toFixed(4)} mmol product per mmol substrate (product floor ${fnum(res.floor, 4)} over ${esc(subEx)} uptake ${y.uptake.toFixed(3)} ${UNIT})` +
            (y.cmol != null ? ` · ${y.cmol.toFixed(4)} C-mol per C-mol (from the recorded formulas)` : ` · ${esc(y.note || '')}`)}.
        Flux ratios from one LP solution; model predictions, not measurements.</p>`;
    } else if (res.verdict === 'not-coupled') {
      numbers = `<div class="an-chips">
        ${chipHTML('na', `max growth ${fnum(res.mu, 4)} h<sup>-1</sup>`)}
        ${res.floor != null ? chipHTML('bad', `floor ${fnum(res.floor, 4)} ${UNIT}`) : ''}
        ${res.ratio != null ? chipHTML('na', `decoupled-growth ratio ${fnum(res.ratio, 3)} (coupling needs < ${res.frac})`) : ''}
      </div>`;
    }
    const mapBtn = res.koSet.length && ctx.mapAvailable()
      ? `<button class="btn small" type="button" id="an-couple-map">Show knockouts on the 3D map</button>` : '';
    const applyBtn = res.koSet.length
      ? `<button class="btn small" type="button" id="an-couple-apply">Apply knockout set to the session</button>` : '';
    out.innerHTML = `
      <div class="an-chips" style="margin-top:var(--s3)">${verdictChip}</div>
      <p class="status">${esc(res.reason)} Search used ${res.solves} LPs in ${secs} s.</p>
      ${numbers}
      ${koTable}
      <div class="cardactions">${applyBtn} ${mapBtn}</div>
      ${iterLines.length ? `<details style="margin-top:var(--s3)"><summary style="cursor:pointer">Search log (${iterLines.length} iteration${iterLines.length === 1 ? '' : 's'})</summary>
        <ul class="an-log">${iterLines.map(l => `<li>${esc(l)}</li>`).join('')}</ul></details>` : ''}`;
    const apply = out.querySelector('#an-couple-apply');
    if (apply) apply.addEventListener('click', () => {
      for (const k of res.koSet) setEdit(st.acc, k.id, 0, 0);
      apply.disabled = true;
      apply.textContent = `Applied ${res.koSet.length} knockout${res.koSet.length === 1 ? '' : 's'} (see chips above)`;
    });
    const mapB = out.querySelector('#an-couple-map');
    if (mapB) mapB.addEventListener('click', () => {
      const r = ctx.showKOsOnMap(res.koSet.map(k => k.id));
      mapB.insertAdjacentHTML('afterend', `<span class="status"> ${r
        ? `Marked ${r.drawn} of ${res.koSet.length} knockout reaction${res.koSet.length === 1 ? '' : 's'} on the union map${r.missing.length ? ` (${r.missing.length} not in the union graph)` : ''}.`
        : 'The 3D map is unavailable in this session.'}</span>`);
    });
  }

  // ------------------------------------------------------------- envelope ----
  $('#an-env-cancel').addEventListener('click', () => { st.cancelFlag = true; });
  $('#an-env-run').addEventListener('click', async () => {
    if (!ready() || !st.prod || st.busy) return;
    startRun('envelope');
    const prog = $('#an-env-prog');
    const btnCancel = $('#an-env-cancel');
    btnCancel.hidden = false;
    const n = clampInt($('#an-envn').value, 8, 60, 24);
    try {
      const { S, target } = await newSession();
      if (!target) {
        $('#an-env-out').innerHTML = `<p class="status error">${esc(st.acc)} does not contain the chosen product.</p>`;
        return;
      }
      prog.textContent = `Solving ${2 * n + 1} LPs…`;
      const wt = await productionEnvelope(S, n, {
        shouldStop: () => st.cancelFlag,
        onProgress: (d, t) => { prog.textContent = `Envelope point ${d} of ${t} (2 LPs each)…`; },
      });
      if (!wt.ok) {
        $('#an-env-out').innerHTML = `<p class="status error">Envelope not computed: the growth LP was ${esc(wt.statusText)} on this medium.</p>`;
        return;
      }
      let ko = null;
      const koSet = (st.search && st.search.res.koSet.length && st.search.stamp === staleStamp())
        ? st.search.res.koSet.map(k => k.id) : null;
      const koApplied = new Set(sessionKOs());
      const overlayIds = koSet && koSet.some(id => !koApplied.has(id)) ? koSet : null;
      if (overlayIds && !st.cancelFlag) {
        prog.textContent = `Reference envelope done; solving the knockout envelope (${2 * n + 1} LPs)…`;
        const { S: S2 } = await newSession(overlayIds);
        ko = await productionEnvelope(S2, n, {
          shouldStop: () => st.cancelFlag,
          onProgress: (d, t) => { prog.textContent = `Knockout envelope point ${d} of ${t}…`; },
        });
        if (!ko.ok) ko = null;
      }
      st.env = { wt, ko, koIds: ko ? overlayIds : null, n, stamp: staleStamp() };
      prog.textContent = '';
      renderEnvelope();
    } catch (e) {
      prog.textContent = '';
      $('#an-env-out').innerHTML = `<p class="status error">Envelope failed: ${esc(e.message)}.</p>`;
    } finally {
      btnCancel.hidden = true;
      endRun();
    }
  });

  function renderEnvelope() {
    const host = $('#an-env-out');
    const { wt, ko, koIds, n } = st.env;
    const pts = wt.points.filter(p => p.min != null && p.max != null);
    const failed = wt.points.length - pts.length;
    if (pts.length < 2) {
      host.innerHTML = `<p class="status error">Only ${pts.length} of ${wt.points.length} envelope points solved; nothing to draw.</p>`;
      return;
    }
    const koPts = ko ? ko.points.filter(p => p.min != null && p.max != null) : null;
    const allY = pts.flatMap(p => [p.min, p.max]).concat(koPts ? koPts.flatMap(p => [p.min, p.max]) : []);
    const allX = pts.map(p => p.mu).concat(koPts ? koPts.map(p => p.mu) : []);
    const svg = envelopeSVG(pts, koPts, Math.max(...allX), Math.min(0, Math.min(...allY)), Math.max(...allY, FLUX_TOL));
    const prodName = ctx.metName(st.prod) || st.prod;
    host.innerHTML = `
      ${svg}
      <p class="status">Production envelope of ${esc(prodName)} for ${esc(st.acc)} on ${esc(st.mediumLabel)}:
      ${pts.length} of ${wt.points.length} biomass levels solved${failed ? ` (${failed} failed and are not drawn)` : ''},
      ${wt.solves} LPs. The reference max growth is marked; the band spans the minimum to maximum product flux at each biomass level.
      ${koPts ? `The knockout envelope (${koIds.map(esc).join(', ')}) is the accent band; its minimum staying above zero at high biomass is the growth coupling.` : ''}
      ${!koPts && st.search && st.search.res.koSet.length && sessionKOs().length ? 'The searched knockout set is already applied to the session, so the single band above includes it.' : ''}</p>
      <div class="cardactions"><button class="btn small" type="button" id="an-env-csv">Export points CSV</button></div>`;
    host.querySelector('#an-env-csv').addEventListener('click', () => {
      const lines = [
        `# production envelope: GEM ${st.acc} on ${st.mediumLabel}; product ${st.prod}; ${pts.length} of ${wt.points.length} points solved; units mmol/gDW/h (biomass 1/h)`,
        'series,biomass,product_min,product_max',
        ...wt.points.map(p => ['reference', p.mu, p.min ?? '', p.max ?? ''].map(csvEscape).join(',')),
        ...(ko ? ko.points.map(p => [`knockout:${koIds.join('+')}`, p.mu, p.min ?? '', p.max ?? ''].map(csvEscape).join(',')) : []),
      ];
      downloadBlob(lines.join('\n'), `envelope_${st.acc}_${st.prod}.csv`, 'text/csv');
    });
  }

  function envelopeSVG(pts, koPts, xMax, yMin, yMax) {
    const W = 720, H = 400, mL = 64, mR = 16, mT = 16, mB = 48;
    const iw = W - mL - mR, ih = H - mT - mB;
    if (xMax <= 0) xMax = 1;
    const ySpan = (yMax - yMin) || 1;
    const X = (v) => mL + v / xMax * iw;
    const Y = (v) => mT + (yMax - v) / ySpan * ih;
    const ticksX = niceTicks(0, xMax, 6);
    const ticksY = niceTicks(yMin, yMax, 6);
    const line = (arr, key) => arr.map((p, i) => `${i ? 'L' : 'M'}${X(p.mu).toFixed(1)},${Y(p[key]).toFixed(1)}`).join(' ');
    const band = (arr) => line(arr, 'max') + ' ' + [...arr].reverse().map(p => `L${X(p.mu).toFixed(1)},${Y(p.min).toFixed(1)}`).join(' ') + ' Z';
    const muMark = pts[pts.length - 1].mu;
    return `
    <svg class="an-envsvg" viewBox="0 0 ${W} ${H}" role="img"
      aria-label="Production envelope: product flux versus biomass flux, minimum and maximum product at each biomass level${koPts ? ', reference and knockout series' : ''}.">
      <rect x="${mL}" y="${mT}" width="${iw}" height="${ih}" fill="none" stroke="var(--line)"/>
      ${ticksY.map(t => `<line x1="${mL}" x2="${W - mR}" y1="${Y(t)}" y2="${Y(t)}" stroke="var(--line)" stroke-width="0.5"/>
        <text x="${mL - 6}" y="${Y(t) + 4}" text-anchor="end" class="an-envtick">${tickLabel(t)}</text>`).join('')}
      ${ticksX.map(t => `<line y1="${mT}" y2="${H - mB}" x1="${X(t)}" x2="${X(t)}" stroke="var(--line)" stroke-width="0.5"/>
        <text x="${X(t)}" y="${H - mB + 16}" text-anchor="middle" class="an-envtick">${tickLabel(t)}</text>`).join('')}
      ${yMin < 0 ? `<line x1="${mL}" x2="${W - mR}" y1="${Y(0)}" y2="${Y(0)}" stroke="var(--line-strong)"/>` : ''}
      <path d="${band(pts)}" fill="${koPts ? 'var(--surface-2)' : 'var(--accent-wash)'}" stroke="none"/>
      <path d="${line(pts, 'max')}" fill="none" stroke="${koPts ? 'var(--ink-2)' : 'var(--accent-ink)'}" stroke-width="1.6"/>
      <path d="${line(pts, 'min')}" fill="none" stroke="${koPts ? 'var(--ink-2)' : 'var(--accent-ink)'}" stroke-width="1.6" stroke-dasharray="${koPts ? '' : '5 3'}"/>
      ${koPts ? `
        <path d="${band(koPts)}" fill="var(--accent-wash)" fill-opacity="0.85" stroke="none"/>
        <path d="${line(koPts, 'max')}" fill="none" stroke="var(--accent-ink)" stroke-width="1.8"/>
        <path d="${line(koPts, 'min')}" fill="none" stroke="var(--accent-ink)" stroke-width="1.8" stroke-dasharray="5 3"/>` : ''}
      <line x1="${X(muMark)}" x2="${X(muMark)}" y1="${mT}" y2="${H - mB}" stroke="var(--line-strong)" stroke-dasharray="2 3"/>
      <text x="${Math.min(X(muMark), W - mR - 4)}" y="${mT + 14}" text-anchor="end" class="an-envtick">max growth ${fnum(muMark, 3)}</text>
      <text x="${mL + iw / 2}" y="${H - 8}" text-anchor="middle" class="an-envlabel">biomass flux (h⁻¹)</text>
      <text transform="translate(14 ${mT + ih / 2}) rotate(-90)" text-anchor="middle" class="an-envlabel">product flux (mmol gDW⁻¹ h⁻¹)</text>
      ${koPts ? `
        <g class="an-envtick">
          <rect x="${mL + 10}" y="${mT + 8}" width="14" height="10" fill="var(--surface-2)" stroke="var(--ink-2)"/>
          <text x="${mL + 30}" y="${mT + 17}">reference</text>
          <rect x="${mL + 10}" y="${mT + 24}" width="14" height="10" fill="var(--accent-wash)" stroke="var(--accent-ink)"/>
          <text x="${mL + 30}" y="${mT + 33}">with knockouts</text>
        </g>` : ''}
    </svg>`;
  }

  function tickLabel(t) {
    const a = Math.abs(t);
    if (a >= 100 || a === 0) return String(Math.round(t));
    if (a >= 1) return t.toFixed(1).replace(/\.0$/, '');
    return t.toPrecision(2);
  }
  function niceTicks(lo, hi, n) {
    const span = hi - lo || 1;
    const step0 = span / n;
    const mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const step = [1, 2, 5, 10].map(m => m * mag).find(s => span / s <= n) || 10 * mag;
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(v);
    return out;
  }

  // -------------------------------------------------------------- exports ----
  renderSweepTable();
  return {
    setActive(v) {
      if (v && st.gem) prefillEndpoints();
    },
  };
}
