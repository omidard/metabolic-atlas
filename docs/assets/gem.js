// Model stage: the chosen GEM as a dashboard (composition, balance, biomass,
// maintenance, growth anchors) over the browsable tables (reactions with
// editable bounds, metabolites, genes with sequences). One GEM at a time,
// lazily; the choice is written to the session context and carries forward.

import { loadGem, loadSeqs, getEdit, setEdit, clearEdit, clearAllEdits, editCount, fmt, downloadBlob, csvEscape, DATA_RELEASE } from './data.js';
import { getContext, setContext, onContext } from './context.js';
import { GROUP_COLORS, NEUTRAL_BAR, chartBlock, hBars, stackBar, donut, fractionBar, intervals, statCard, fitWidth } from './charts.js';

const esc = (s) => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

let root, statusEl, gemSelect;
let gem = null;               // loaded GEM json
let acc = null;
let seqs = null;              // loaded sequences for acc
let rxnFilter = '', rxnPage = 0, rxnPageSize = 25;
let metFilter = '', metPage = 0;
let geneFilter = '';
let pendingRxnFilter = null;
let onAccentChange = null;

export function initGemView(container, indexData, opts = {}) {
  root = container;
  onAccentChange = opts.onAccentChange || null;
  root.innerHTML = `
    <h1>Which model is this, and what is in it?</h1>
    <p class="sub">35 genome-scale models across 4 species. Choosing a model here carries it into
    Simulate and Engineer. The dashboard summarises composition, mass balance, biomass and growth
    anchors; the tables below it hold every reaction (bounds editable), metabolite and gene
    sequence. Session bound edits apply to every solve and are written into the SBML and COBRA
    JSON exports.</p>
    <div class="gem-toolbar">
      <div class="field">
        <label for="gem-select">Model (35 GEMs)</label>
        <select id="gem-select"><option value="">Choose a GEM…</option></select>
      </div>
      <span id="gem-status" class="status" role="status" aria-live="polite"></span>
    </div>
    <div class="card empty-state" id="gem-empty">
      <strong>No model chosen yet.</strong>
      <p style="margin:var(--s2) 0 0">Choose one of the 35 GEMs above, or pick a strain from a
      Discover result. The dashboard renders from the model file (about 1 MB, loaded once).</p>
    </div>
    <div id="gem-body" hidden></div>`;
  statusEl = root.querySelector('#gem-status');
  gemSelect = root.querySelector('#gem-select');

  const bySpecies = new Map();
  for (const g of indexData.gems) {
    if (!bySpecies.has(g.species)) bySpecies.set(g.species, []);
    bySpecies.get(g.species).push(g);
  }
  for (const [sp, gems] of bySpecies) {
    const og = document.createElement('optgroup');
    og.label = `${sp} (${gems.length})`;
    for (const g of gems) {
      const o = document.createElement('option');
      o.value = g.acc;
      o.textContent = `${g.acc} · ${fmt.format(g.reactions)} rxns · ${fmt.format(g.genes)} genes`;
      og.appendChild(o);
    }
    gemSelect.appendChild(og);
  }
  gemSelect.addEventListener('change', () => selectGem(gemSelect.value));

  // context sync: a GEM chosen in Simulate or Engineer appears here too
  onContext((c, changed) => {
    if (changed.includes('gem') && (c.gem || '') !== (acc || '') && (c.gem || '') !== gemSelect.value) {
      gemSelect.value = c.gem || '';
      selectGem(c.gem || '');
    }
  });
  const c0 = getContext();
  if (c0.gem) { gemSelect.value = c0.gem; selectGem(c0.gem); }
}

export function openReactionInBrowser(rxnId) {
  pendingRxnFilter = rxnId;
  if (gem) {
    rxnFilter = rxnId; rxnPage = 0;
    const inp = root.querySelector('#rxn-search');
    if (inp) { inp.value = rxnId; renderRxnTable(); }
    pendingRxnFilter = null;
  }
}

async function selectGem(newAcc) {
  const body = root.querySelector('#gem-body');
  const empty = root.querySelector('#gem-empty');
  if (!newAcc) { body.hidden = true; if (empty) empty.hidden = false; gem = null; acc = null; setContext({ gem: null }); return; }
  try {
    statusEl.textContent = `Loading GEM ${newAcc} (about 1 MB)…`;
    statusEl.classList.remove('error');
    const g = await loadGem(newAcc);
    gem = g; acc = newAcc; seqs = null;
    rxnFilter = pendingRxnFilter || ''; rxnPage = 0; metFilter = ''; metPage = 0; geneFilter = '';
    pendingRxnFilter = null;
    statusEl.textContent = `${g.species} · ${newAcc}`;
    if (onAccentChange) onAccentChange(g.species);
    setContext({ gem: newAcc });
    renderGem();
    body.hidden = false;
    if (empty) empty.hidden = true;
  } catch (e) {
    statusEl.textContent = `Could not load GEM ${newAcc}: ${e.message}. Check the connection and choose the model again.`;
    statusEl.classList.add('error');
    body.hidden = true;
    if (empty) empty.hidden = false;
  }
}

function renderGem() {
  const s = gem.stats;
  const body = root.querySelector('#gem-body');
  const nBiomass = Object.keys(gem.biomass || {}).length;

  body.innerHTML = `
    <h2>${esc(gem.species)} <span class="mono">${esc(acc)}</span></h2>
    <div class="cardactions" style="margin:0 0 var(--s3)">
      <button class="btn small" id="gem-export-sbml" type="button">Export SBML (FBC v2)</button>
      <button class="btn small" id="gem-export-cobra" type="button">Export COBRA JSON</button>
      <span class="status" id="gem-export-status" role="status" aria-live="polite">Whole-model exports; session bound edits are written into the file.</span>
    </div>
    ${dashboardHTML()}

    <details style="margin-top:var(--s4)">
      <summary style="cursor:pointer;font-weight:600">Biomass reaction (${esc(s.biomass_id)}, ${fmt.format(nBiomass)} components)</summary>
      <p class="biomass-formula" id="biomass-formula"></p>
    </details>

    <div class="tabs" role="tablist" aria-label="GEM sections">
      <button class="tab" role="tab" id="tab-rxns" aria-selected="true" aria-controls="panel-rxns">Reactions (${fmt.format(gem.reactions.length)})</button>
      <button class="tab" role="tab" id="tab-mets" aria-selected="false" aria-controls="panel-mets">Metabolites (${fmt.format(gem.metabolites.length)})</button>
      <button class="tab" role="tab" id="tab-genes" aria-selected="false" aria-controls="panel-genes">Genes (${fmt.format(gem.genes.length)})</button>
    </div>
    <div id="panel-rxns" class="tabpanel" role="tabpanel" aria-labelledby="tab-rxns"></div>
    <div id="panel-mets" class="tabpanel" role="tabpanel" aria-labelledby="tab-mets" hidden></div>
    <div id="panel-genes" class="tabpanel" role="tabpanel" aria-labelledby="tab-genes" hidden></div>`;

  body.querySelector('#gem-export-sbml').addEventListener('click', () => exportModel('sbml'));
  body.querySelector('#gem-export-cobra').addEventListener('click', () => exportModel('cobra'));

  renderBiomass();
  renderRxnPanel();
  renderMetPanel();
  renderGenePanel();

  const tabs = body.querySelectorAll('.tab');
  tabs.forEach(t => t.addEventListener('click', () => {
    tabs.forEach(x => x.setAttribute('aria-selected', String(x === t)));
    body.querySelector('#panel-rxns').hidden = t.id !== 'tab-rxns';
    body.querySelector('#panel-mets').hidden = t.id !== 'tab-mets';
    body.querySelector('#panel-genes').hidden = t.id !== 'tab-genes';
  }));
}

function stat(k, v, d) {
  return `<div class="stat"><div class="k">${k}</div><div class="v">${v}</div><div class="d">${d}</div></div>`;
}

// ---------------- dashboard (all values from the loaded model file) ----------------
const COMP_LABEL = { c: 'cytosol', p: 'periplasm', e: 'extracellular' };
const COMP_COLORS = { cytosol: '#3A62B0', periplasm: '#2A8A80', extracellular: '#A65D1E', 'cross-compartment': '#8A4BA8', exchange: '#98948C' };

function normComp(cRaw) {
  // compartment ids are c/p/e, occasionally prefixed (C_c); take the suffix
  const c = String(cRaw || '');
  const t = c.includes('_') ? c.split('_').pop() : c;
  return COMP_LABEL[t] || (c || 'unrecorded');
}

function dashboardHTML() {
  const s = gem.stats;
  const R = gem.reactions.length;

  // headline cards: this model versus the mean of all registered GEMs
  const allGems = indexCache ? indexCache.gems : [];
  const mean = (key) => allGems.length ? allGems.reduce((a, g) => a + g[key], 0) / allGems.length : null;
  const meanLabel = `${allGems.length}-GEM mean`;
  const speciesAccent = (indexCache && (indexCache.gems.find(g => g.acc === acc) || {}).accent) || 'var(--accent)';

  // 1 · reactions by pathway group
  const byGroup = new Map();
  for (const r of gem.reactions) {
    const g = r.group || 'Ungrouped';
    byGroup.set(g, (byGroup.get(g) || 0) + 1);
  }
  const groupRows = [...byGroup.entries()]
    .sort((a, b) => b[1] - a[1])
    .map(([g, n]) => ({ label: g, value: n, color: g === 'Ungrouped' ? '#C9C5BD' : (GROUP_COLORS[g] || NEUTRAL_BAR) }));

  // 2 · reactions by compartment (internal single-compartment, cross-compartment, exchange)
  const metComp = new Map(gem.metabolites.map(m => [m.id, normComp(m.comp)]));
  const compCounts = { cytosol: 0, periplasm: 0, extracellular: 0, 'cross-compartment': 0, exchange: 0 };
  for (const r of gem.reactions) {
    if (r.ex) { compCounts.exchange++; continue; }
    const comps = new Set(Object.keys(r.stoich).map(m => metComp.get(m) || 'unrecorded'));
    if (comps.size === 1) {
      const c = [...comps][0];
      if (c in compCounts) compCounts[c]++; else compCounts['cross-compartment']++;
    } else {
      compCounts['cross-compartment']++;
    }
  }
  const compSegs = Object.entries(compCounts).filter(([, n]) => n > 0)
    .map(([label, n]) => ({ label, value: n, color: COMP_COLORS[label] }));

  // 3 · metabolites by compartment
  const metCompCounts = new Map();
  for (const m of gem.metabolites) {
    const c = normComp(m.comp);
    metCompCounts.set(c, (metCompCounts.get(c) || 0) + 1);
  }
  const metSegs = [...metCompCounts.entries()].sort((a, b) => b[1] - a[1])
    .map(([label, n]) => ({ label, value: n, color: COMP_COLORS[label] || NEUTRAL_BAR }));

  // 5 · exchange / transport / internal
  const nTransport = gem.reactions.filter(r => r.transport && !r.ex).length;
  const roleSegs = [
    { label: 'internal', value: R - s.exchanges - nTransport, color: speciesAccent },
    { label: 'transport', value: nTransport, color: '#2E7A9E' },
    { label: 'exchange', value: s.exchanges, color: '#98948C' },
  ];

  // 7 · biomass precursors (consumed side, maintenance ATP terms set aside)
  const GAM_TERMS = new Set(['atp_c', 'h2o_c', 'adp_c', 'pi_c', 'h_c']);
  const consumed = Object.entries(gem.biomass || {}).filter(([, c]) => c < 0);
  const precursors = consumed.filter(([m]) => !GAM_TERMS.has(m))
    .sort((a, b) => Math.abs(b[1]) - Math.abs(a[1]));
  const topPre = precursors.slice(0, 12).map(([m, c]) => ({
    label: m, value: Math.abs(c), color: speciesAccent,
  }));

  // 9 · growth anchors
  const growthRows = Object.entries(gem.growth || {}).map(([axis, r]) => ({
    label: axis, lo: r?.[0] ?? null, hi: r?.[1] ?? null,
  }));

  const w = fitWidth(root, 380);
  return `
    <div class="dash-head">
      ${statCard('Genes', s.genes, 'in model', { mean: mean('genes'), meanLabel, color: speciesAccent })}
      ${statCard('Reactions', s.reactions, 'in model', { mean: mean('reactions'), meanLabel, color: speciesAccent })}
      ${statCard('Metabolites', s.metabolites, 'in model', { mean: mean('metabolites'), meanLabel, color: speciesAccent })}
      ${statCard('Exchanges', s.exchanges, `of ${fmt.format(R)} reactions`)}
      ${statCard('Transporters', s.transporters, `of ${fmt.format(R)} reactions`)}
    </div>

    <div class="chart-grid">
      ${chartBlock(`Reaction composition by pathway group (${fmt.format(R)} reactions)`,
        hBars(groupRows, { total: R }),
        'Groups follow the union-map assignment; reactions outside the 13 named groups are Ungrouped.')}

      ${chartBlock(`Reactions by compartment (${fmt.format(R)} reactions)`,
        stackBar(compSegs, { total: R }),
        'Internal reactions are placed by the compartment of their metabolites; a reaction spanning two compartments counts as cross-compartment.')}

      ${chartBlock(`Metabolites by compartment (${fmt.format(gem.metabolites.length)} metabolites)`,
        donut(metSegs, { centerTop: fmt.format(gem.metabolites.length), centerBottom: 'metabolites' }))}

      ${chartBlock(`Exchange, transport and internal reactions (${fmt.format(R)} reactions)`,
        stackBar(roleSegs, { total: R }))}

      ${chartBlock('Mass balance',
        fractionBar(s.mass_balanced, s.balanceable, { color: speciesAccent, label: 'mass balanced' }),
        `Balanced of the ${fmt.format(s.balanceable)} balanceable reactions (${fmt.format(R - s.balanceable)} of ${fmt.format(R)} carry no full formula set and cannot be checked).`)}

      ${chartBlock('Reactions without a gene rule (gap-filled or orphan)',
        s.gapfilled_orphan == null
          ? '<p class="status">not computed for this model</p>'
          : fractionBar(s.gapfilled_orphan, R, { color: '#A65D1E', label: 'without a gene rule' }),
        s.gapfilled_orphan == null ? '' : 'Exchange, biomass and demand reactions are excluded from the count.')}

      ${chartBlock(`Biomass precursor coefficients (top ${topPre.length} of ${fmt.format(precursors.length)} consumed components)`,
        hBars(topPre, { showPct: false, valueFmt: (v) => String(Math.round(v * 1e4) / 1e4) }),
        `|coefficient| in mmol gDW<sup>-1</sup>. The maintenance terms (ATP, H<sub>2</sub>O, ADP, Pi, H; GAM ${s.gam ?? 'not computed'}) are set aside; the full ${fmt.format(Object.keys(gem.biomass || {}).length)}-component reaction is below.`)}

      ${chartBlock('Maintenance energy and biomass objective',
        `<div class="dash-head" style="margin-top:0">
          ${stat('GAM', s.gam ?? 'not computed', 'mmol ATP gDW<sup>-1</sup>, growth-associated')}
          ${stat('NGAM', s.ngam ?? 'not computed', `mmol ATP gDW<sup>-1</sup> h<sup>-1</sup> · ${esc(s.ngam_source || 'source not recorded')}`)}
          ${stat('Biomass', `<span class="mono" style="font-size:var(--t-s)">${esc(s.biomass_id)}</span>`, `${fmt.format(Object.keys(gem.biomass || {}).length)} components`)}
        </div>`)}

      ${chartBlock(`Growth anchors per medium axis (${growthRows.length} ${growthRows.length === 1 ? 'axis' : 'axes'}, model-derived)`,
        growthRows.length
          ? intervals(growthRows, { width: Math.min(w, 380), unit: 'h^-1', color: speciesAccent })
          : '<p class="status">No growth anchors recorded for this model.</p>',
        'Anchor interval (mu min to mu max, h<sup>-1</sup>) recorded in the model release; a reference interval from the model chain, not an experimental measurement.')}
    </div>`;
}

// Whole-model export: SBML (level 3 version 1, fbc v2) or COBRA JSON, both
// with the session's in-memory bound edits applied at click time. The builder
// loads lazily on first use.
async function exportModel(kind) {
  const el = root.querySelector('#gem-export-status');
  try {
    el.textContent = 'Building the export…';
    const mod = await import('./sbml.js');
    const nEdits = editCount(acc);
    if (kind === 'sbml') {
      const { xml, nGprSkipped } = mod.buildSBML(gem, acc, { release: DATA_RELEASE });
      downloadBlob(xml, `${acc}.xml`, 'application/xml');
      el.textContent = `${acc}.xml written: ${fmt.format(gem.reactions.length)} reactions, ` +
        `${fmt.format(gem.metabolites.length)} metabolites, ${fmt.format(gem.genes.length)} gene products; ` +
        `${nEdits} of ${fmt.format(gem.reactions.length)} reaction bounds carry session edits` +
        (nGprSkipped ? `; ${nGprSkipped} gene rules could not be parsed and are omitted (kept in the JSON export)` : '') + '.';
    } else {
      const { obj } = mod.buildCobraJSON(gem, acc, { release: DATA_RELEASE });
      downloadBlob(JSON.stringify(obj), `${acc}.cobra.json`, 'application/json');
      el.textContent = `${acc}.cobra.json written: ${fmt.format(gem.reactions.length)} reactions; ` +
        `${nEdits} of ${fmt.format(gem.reactions.length)} reaction bounds carry session edits.`;
    }
  } catch (e) {
    el.textContent = `Export failed: ${e.message}. Choose the button again to retry.`;
  }
}

function renderBiomass() {
  const el = root.querySelector('#biomass-formula');
  const terms = Object.entries(gem.biomass || {});
  const lhs = terms.filter(([, c]) => c < 0).map(([m, c]) => `${fmtCoef(-c)} ${m}`).join(' + ');
  const rhs = terms.filter(([, c]) => c > 0).map(([m, c]) => `${fmtCoef(c)} ${m}`).join(' + ');
  el.textContent = `${lhs} → ${rhs || 'biomass'}`;
}

const fmtCoef = (c) => (c === 1 ? '' : String(Math.round(c * 1e6) / 1e6)).trim() || '1';

// ---------------- reactions ----------------
function renderRxnPanel() {
  const p = root.querySelector('#panel-rxns');
  p.innerHTML = `
    <div class="tablebar">
      <div class="field" style="flex:1 1 220px">
        <label for="rxn-search">Filter reactions (id, name, subsystem, GPR, EC)</label>
        <input type="search" id="rxn-search" value="${esc(rxnFilter)}" placeholder="e.g. PGK or Glycolysis">
      </div>
      <div class="field">
        <label for="rxn-pagesize">Rows per page</label>
        <select id="rxn-pagesize"><option>25</option><option>50</option><option>100</option></select>
      </div>
      <button class="btn small" id="rxn-export-csv" type="button">Export CSV</button>
      <button class="btn small" id="rxn-export-json" type="button">Export JSON</button>
      <button class="btn small" id="rxn-reset-all" type="button">Reset all edits</button>
    </div>
    <p class="count" id="rxn-count" role="status" aria-live="polite"></p>
    <div class="tablewrap"><table class="data" id="rxn-table">
      <thead><tr>
        <th scope="col">Id</th><th scope="col">Name</th><th scope="col">Subsystem</th>
        <th scope="col">Equation</th><th scope="col">lb</th><th scope="col">ub</th>
        <th scope="col">GPR</th><th scope="col">EC / cross-refs</th><th scope="col"></th>
      </tr></thead><tbody></tbody>
    </table></div>
    <div class="pager" style="margin-top:var(--s3)">
      <button class="btn small" id="rxn-prev" type="button">Previous</button>
      <span id="rxn-pageinfo"></span>
      <button class="btn small" id="rxn-next" type="button">Next</button>
    </div>`;
  p.querySelector('#rxn-search').addEventListener('input', (e) => { rxnFilter = e.target.value.trim(); rxnPage = 0; renderRxnTable(); });
  p.querySelector('#rxn-pagesize').addEventListener('change', (e) => { rxnPageSize = +e.target.value; rxnPage = 0; renderRxnTable(); });
  p.querySelector('#rxn-prev').addEventListener('click', () => { if (rxnPage > 0) { rxnPage--; renderRxnTable(); } });
  p.querySelector('#rxn-next').addEventListener('click', () => { rxnPage++; renderRxnTable(); });
  p.querySelector('#rxn-export-csv').addEventListener('click', () => exportRxns('csv'));
  p.querySelector('#rxn-export-json').addEventListener('click', () => exportRxns('json'));
  p.querySelector('#rxn-reset-all').addEventListener('click', () => { clearAllEdits(acc); renderRxnTable(); });
  renderRxnTable();
}

function filteredRxns() {
  if (!rxnFilter) return gem.reactions;
  const q = rxnFilter.toLowerCase();
  return gem.reactions.filter(r =>
    r.id.toLowerCase().includes(q) || (r.name || '').toLowerCase().includes(q) ||
    (r.subsystem || '').toLowerCase().includes(q) || (r.gpr || '').toLowerCase().includes(q) ||
    ((r.xr && r.xr.ec) || []).some(e => e.includes(q)));
}

function effBounds(r) {
  const e = getEdit(acc, r.id);
  return e ? { lb: e.lb, ub: e.ub, edited: true } : { lb: r.lb, ub: r.ub, edited: false };
}

function equation(r, b) {
  const lhs = [], rhs = [];
  for (const [m, c] of Object.entries(r.stoich)) {
    (c < 0 ? lhs : rhs).push(`${fmtCoef(Math.abs(c))} ${m}`.trim());
  }
  const arrow = (b.lb < 0 && b.ub > 0) ? '⇌' : (b.ub <= 0 && b.lb < 0) ? '←' : '→';
  return `${lhs.join(' + ')} ${arrow} ${rhs.join(' + ')}`;
}

function xrLinks(r) {
  const x = r.xr || {};
  const out = [];
  for (const ec of x.ec || []) out.push(`<span class="mono">EC ${esc(ec)}</span>`);
  for (const id of x.kegg || []) out.push(`<a href="https://www.kegg.jp/entry/${esc(id)}" target="_blank" rel="noopener">KEGG ${esc(id)}</a>`);
  for (const id of x.bigg || []) out.push(`<a href="https://bigg.ucsd.edu/universal/reactions/${esc(id)}" target="_blank" rel="noopener">BiGG</a>`);
  for (const id of x.rhea || []) out.push(`<a href="https://www.rhea-db.org/rhea/${esc(id)}" target="_blank" rel="noopener">RHEA ${esc(id)}</a>`);
  for (const id of x.metanetx || []) out.push(`<a href="https://www.metanetx.org/equa_info/${esc(id)}" target="_blank" rel="noopener">MetaNetX</a>`);
  return out.join(' ') || '<span class="status">none</span>';
}

function renderRxnTable() {
  const list = filteredRxns();
  const pages = Math.max(1, Math.ceil(list.length / rxnPageSize));
  rxnPage = Math.min(rxnPage, pages - 1);
  const slice = list.slice(rxnPage * rxnPageSize, (rxnPage + 1) * rxnPageSize);

  const nEdits = editCount(acc);
  root.querySelector('#rxn-count').textContent =
    (rxnFilter
      ? `Showing ${fmt.format(slice.length)} of ${fmt.format(list.length)} matching reactions (filtered from ${fmt.format(gem.reactions.length)})`
      : `Showing ${fmt.format(slice.length)} of ${fmt.format(list.length)} reactions`) +
    (nEdits ? ` · ${nEdits} bound edit${nEdits > 1 ? 's' : ''} in memory` : '');
  root.querySelector('#rxn-pageinfo').textContent = `Page ${rxnPage + 1} of ${pages}`;
  root.querySelector('#rxn-prev').disabled = rxnPage === 0;
  root.querySelector('#rxn-next').disabled = rxnPage >= pages - 1;

  const tbody = root.querySelector('#rxn-table tbody');
  tbody.innerHTML = slice.map(r => {
    const b = effBounds(r);
    return `<tr data-rid="${esc(r.id)}" class="${b.edited ? 'edited' : ''}">
      <td class="mono">${esc(r.id)}</td>
      <td>${esc(r.name || '')}</td>
      <td>${esc(r.subsystem || '')}</td>
      <td class="eq">${esc(equation(r, b))}</td>
      <td><input class="bound" type="number" step="any" value="${b.lb}" aria-label="Lower bound of ${esc(r.id)}" data-kind="lb"></td>
      <td><input class="bound" type="number" step="any" value="${b.ub}" aria-label="Upper bound of ${esc(r.id)}" data-kind="ub"></td>
      <td class="mono" style="max-width:220px;overflow-wrap:anywhere">${esc(r.gpr || '')}</td>
      <td class="xr">${xrLinks(r)}</td>
      <td>${b.edited ? `<button class="btn small" type="button" data-reset="${esc(r.id)}">Reset</button>` : ''}</td>
    </tr>`;
  }).join('');

  tbody.querySelectorAll('input.bound').forEach(inp => {
    inp.addEventListener('change', () => {
      const tr = inp.closest('tr');
      const rid = tr.dataset.rid;
      const r = gem.reactions.find(x => x.id === rid);
      const cur = effBounds(r);
      const lb = inp.dataset.kind === 'lb' ? parseFloat(inp.value) : cur.lb;
      const ub = inp.dataset.kind === 'ub' ? parseFloat(inp.value) : cur.ub;
      if (Number.isNaN(lb) || Number.isNaN(ub)) return;
      if (lb === r.lb && ub === r.ub) clearEdit(acc, rid); else setEdit(acc, rid, lb, ub);
      renderRxnTable();
    });
  });
  tbody.querySelectorAll('button[data-reset]').forEach(btn => {
    btn.addEventListener('click', () => { clearEdit(acc, btn.dataset.reset); renderRxnTable(); });
  });
}

function exportRxns(kind) {
  const list = filteredRxns();
  const rows = list.map(r => {
    const b = effBounds(r);
    return {
      id: r.id, name: r.name || '', subsystem: r.subsystem || '', group: r.group || '',
      lb: b.lb, ub: b.ub, edited: b.edited, gpr: r.gpr || '',
      ec: ((r.xr && r.xr.ec) || []).join(';'), equation: equation(r, b),
    };
  });
  const nEdits = rows.filter(r => r.edited).length;
  const stamp = `${acc} reactions: ${rows.length} of ${gem.reactions.length}` + (rxnFilter ? ` (filter: ${rxnFilter})` : '') + `; ${nEdits} edited bounds`;
  if (kind === 'csv') {
    const head = 'id,name,subsystem,group,lb,ub,edited,gpr,ec,equation';
    const csv = [`# ${stamp}`, head, ...rows.map(r => [r.id, r.name, r.subsystem, r.group, r.lb, r.ub, r.edited, r.gpr, r.ec, r.equation].map(csvEscape).join(','))].join('\n');
    downloadBlob(csv, `${acc}_reactions.csv`, 'text/csv');
  } else {
    downloadBlob(JSON.stringify({ note: stamp, reactions: rows }, null, 1), `${acc}_reactions.json`, 'application/json');
  }
}

// ---------------- metabolites ----------------
function metXr(m) {
  const x = m.xr || {};
  const out = [];
  for (const id of x.kegg || []) out.push(`<a href="https://www.kegg.jp/entry/${esc(id)}" target="_blank" rel="noopener">KEGG</a>`);
  for (const id of x.chebi || []) out.push(`<a href="https://www.ebi.ac.uk/chebi/searchId.do?chebiId=${esc(id)}" target="_blank" rel="noopener">ChEBI</a>`);
  for (const id of x.bigg || []) out.push(`<a href="https://bigg.ucsd.edu/universal/metabolites/${esc(id)}" target="_blank" rel="noopener">BiGG</a>`);
  for (const id of x.metanetx || []) out.push(`<a href="https://www.metanetx.org/chem_info/${esc(id)}" target="_blank" rel="noopener">MetaNetX</a>`);
  for (const id of x.seed || []) out.push(`<a href="https://modelseed.org/biochem/compounds/${esc(id)}" target="_blank" rel="noopener">SEED</a>`);
  for (const id of x.hmdb || []) out.push(`<a href="https://hmdb.ca/metabolites/${esc(id)}" target="_blank" rel="noopener">HMDB</a>`);
  for (const id of x.biocyc || []) out.push(`<a href="https://biocyc.org/compound?orgid=META&id=${esc(id)}" target="_blank" rel="noopener">BioCyc</a>`);
  return out.join(' ') || '<span class="status">none</span>';
}

function renderMetPanel() {
  const p = root.querySelector('#panel-mets');
  p.innerHTML = `
    <div class="tablebar">
      <div class="field" style="flex:1 1 220px">
        <label for="met-search">Filter metabolites (id, name, formula)</label>
        <input type="search" id="met-search" placeholder="e.g. glc__D or C6H12O6">
      </div>
    </div>
    <p class="count" id="met-count" role="status" aria-live="polite"></p>
    <div class="tablewrap"><table class="data" id="met-table">
      <thead><tr><th scope="col">Id</th><th scope="col">Name</th><th scope="col">Formula</th>
      <th scope="col">Charge</th><th scope="col">Compartment</th><th scope="col">Currency</th><th scope="col">Cross-refs</th></tr></thead>
      <tbody></tbody></table></div>
    <div class="pager" style="margin-top:var(--s3)">
      <button class="btn small" id="met-prev" type="button">Previous</button>
      <span id="met-pageinfo"></span>
      <button class="btn small" id="met-next" type="button">Next</button>
    </div>`;
  p.querySelector('#met-search').addEventListener('input', (e) => { metFilter = e.target.value.trim().toLowerCase(); metPage = 0; renderMetTable(); });
  p.querySelector('#met-prev').addEventListener('click', () => { if (metPage > 0) { metPage--; renderMetTable(); } });
  p.querySelector('#met-next').addEventListener('click', () => { metPage++; renderMetTable(); });
  renderMetTable();
}

function renderMetTable() {
  const all = gem.metabolites;
  const list = metFilter ? all.filter(m =>
    m.id.toLowerCase().includes(metFilter) || (m.name || '').toLowerCase().includes(metFilter) ||
    (m.formula || '').toLowerCase().includes(metFilter)) : all;
  const size = 25;
  const pages = Math.max(1, Math.ceil(list.length / size));
  metPage = Math.min(metPage, pages - 1);
  const slice = list.slice(metPage * size, (metPage + 1) * size);
  root.querySelector('#met-count').textContent = metFilter
    ? `Showing ${fmt.format(slice.length)} of ${fmt.format(list.length)} matching metabolites (filtered from ${fmt.format(all.length)})`
    : `Showing ${fmt.format(slice.length)} of ${fmt.format(list.length)} metabolites`;
  root.querySelector('#met-pageinfo').textContent = `Page ${metPage + 1} of ${pages}`;
  root.querySelector('#met-prev').disabled = metPage === 0;
  root.querySelector('#met-next').disabled = metPage >= pages - 1;
  root.querySelector('#met-table tbody').innerHTML = slice.map(m => `<tr>
    <td class="mono">${esc(m.id)}</td><td>${esc(m.name || '')}</td>
    <td class="mono">${esc(m.formula || '')}</td><td class="mono">${m.charge ?? ''}</td>
    <td class="mono">${esc(m.comp || '')}</td><td>${m.currency ? 'yes' : 'no'}</td>
    <td class="xr">${metXr(m)}</td></tr>`).join('');
}

// ---------------- genes ----------------
function renderGenePanel() {
  const p = root.querySelector('#panel-genes');
  const nSeqs = (indexEntryFor(acc) || {}).n_seqs;
  p.innerHTML = `
    <p class="count">${fmt.format(gem.genes.length)} genes in the model; nucleotide sequences available for
    ${nSeqs === undefined ? 'not computed' : `${fmt.format(nSeqs)} of ${fmt.format(gem.genes.length)}`}.
    Sequences load on first view (about 1 MB).</p>
    <div class="field" style="max-width:340px">
      <label for="gene-search">Filter genes by id</label>
      <input type="search" id="gene-search" placeholder="e.g. PSEPUT">
    </div>
    <p class="count" id="gene-count" role="status" aria-live="polite"></p>
    <div class="genelist" id="gene-list"></div>
    <div class="seqview" id="seq-view" aria-live="polite"></div>`;
  p.querySelector('#gene-search').addEventListener('input', (e) => { geneFilter = e.target.value.trim().toLowerCase(); renderGeneList(); });
  renderGeneList();
}

let indexCache = null;
export function setIndexData(d) { indexCache = d; }
function indexEntryFor(a) { return indexCache ? indexCache.gems.find(g => g.acc === a) : null; }

function renderGeneList() {
  const list = geneFilter ? gem.genes.filter(g => g.toLowerCase().includes(geneFilter)) : gem.genes;
  const geneMax = 300;
  const slice = list.slice(0, geneMax);
  root.querySelector('#gene-count').textContent =
    slice.length < list.length
      ? `Showing ${fmt.format(slice.length)} of ${fmt.format(list.length)} matching genes; refine the filter to narrow further.`
      : `Showing ${fmt.format(list.length)} of ${fmt.format(gem.genes.length)} genes.`;
  const el = root.querySelector('#gene-list');
  el.innerHTML = slice.map(g => `<button class="btn small" type="button" data-gene="${esc(g)}">${esc(g)}</button>`).join('');
  el.querySelectorAll('button[data-gene]').forEach(b => b.addEventListener('click', () => showSeq(b.dataset.gene)));
}

async function showSeq(geneId) {
  const view = root.querySelector('#seq-view');
  try {
    if (!seqs) {
      view.innerHTML = `<p class="status">Loading sequences for ${esc(acc)} (about 1 MB)…</p>`;
      seqs = await loadSeqs(acc);
    }
    const seq = seqs[geneId];
    if (!seq) {
      view.innerHTML = `<p class="status">No sequence recorded for ${esc(geneId)} in this release.</p>`;
      return;
    }
    const wrapped = seq.replace(/(.{60})/g, '$1\n');
    view.innerHTML = `
      <h3 class="mono">${esc(geneId)} · ${fmt.format(seq.length)} nt</h3>
      <pre>&gt;${esc(acc)}|${esc(geneId)} length=${seq.length}\n${wrapped}</pre>
      <button class="btn small" type="button" id="seq-copy">Copy FASTA</button>`;
    view.querySelector('#seq-copy').addEventListener('click', async () => {
      await navigator.clipboard.writeText(`>${acc}|${geneId} length=${seq.length}\n${wrapped}`);
      view.querySelector('#seq-copy').textContent = 'Copied';
      setTimeout(() => { const b = view.querySelector('#seq-copy'); if (b) b.textContent = 'Copy FASTA'; }, 1500);
    });
  } catch (e) {
    view.innerHTML = `<p class="status error">Could not load sequences for ${esc(acc)}: ${esc(e.message)}. Check the connection and choose the gene again.</p>`;
  }
}
