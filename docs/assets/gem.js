// GEM browser: stats, growth anchors, biomass, reactions (editable bounds),
// metabolites, genes with sequences. One GEM loaded at a time, lazily.

import { loadGem, loadSeqs, getEdit, setEdit, clearEdit, clearAllEdits, editCount, fmt, downloadBlob, csvEscape } from './data.js';

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
    <h1>GEM browser</h1>
    <p class="sub">35 genome-scale models across 4 species. Choose a model to view its
    statistics, growth anchors, biomass formulation, reactions, metabolites and gene
    sequences. Lower and upper bounds are editable in the reactions table; edits are held
    in memory for this session and feed the phase 2 flux analyses.</p>
    <div class="gem-toolbar">
      <div class="field">
        <label for="gem-select">Model (35 GEMs)</label>
        <select id="gem-select"><option value="">Choose a GEM…</option></select>
      </div>
      <span id="gem-status" class="status" role="status" aria-live="polite"></span>
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
  if (!newAcc) { body.hidden = true; gem = null; acc = null; return; }
  try {
    statusEl.textContent = `Loading GEM ${newAcc} (about 1 MB)…`;
    statusEl.classList.remove('error');
    const g = await loadGem(newAcc);
    gem = g; acc = newAcc; seqs = null;
    rxnFilter = pendingRxnFilter || ''; rxnPage = 0; metFilter = ''; metPage = 0; geneFilter = '';
    pendingRxnFilter = null;
    statusEl.textContent = `${g.species} · ${newAcc}`;
    if (onAccentChange) onAccentChange(g.species);
    renderGem();
    body.hidden = false;
  } catch (e) {
    statusEl.textContent = `Could not load GEM ${newAcc}: ${e.message}. Check the connection and choose the model again.`;
    statusEl.classList.add('error');
    body.hidden = true;
  }
}

function renderGem() {
  const s = gem.stats;
  const body = root.querySelector('#gem-body');
  const growthRows = Object.entries(gem.growth || {}).map(([axis, r]) =>
    `<tr><td>${esc(axis)}</td><td class="mono">${r?.[0] ?? 'not computed'}</td><td class="mono">${r?.[1] ?? 'not computed'}</td></tr>`).join('');
  const nBiomass = Object.keys(gem.biomass || {}).length;

  body.innerHTML = `
    <h2>${esc(gem.species)} <span class="mono">${esc(acc)}</span></h2>
    <div class="statgrid">
      ${stat('Genes', fmt.format(s.genes), 'in model')}
      ${stat('Reactions', fmt.format(s.reactions), 'in model')}
      ${stat('Metabolites', fmt.format(s.metabolites), 'in model')}
      ${stat('Exchanges', fmt.format(s.exchanges), `of ${fmt.format(s.reactions)} reactions`)}
      ${stat('Transporters', fmt.format(s.transporters), `of ${fmt.format(s.reactions)} reactions`)}
      ${stat('Mass balanced', fmt.format(s.mass_balanced), `of ${fmt.format(s.reactions)} reactions`)}
      ${stat('Balanceable', fmt.format(s.balanceable), `of ${fmt.format(s.reactions)} reactions`)}
      ${stat('GAM', s.gam ?? 'not computed', 'mmol ATP gDW<sup>-1</sup>')}
      ${stat('NGAM', s.ngam ?? 'not computed', `mmol ATP gDW<sup>-1</sup> h<sup>-1</sup> · ${esc(s.ngam_source || 'source not recorded')}`)}
      ${stat('Biomass', `<span class="mono" style="font-size:var(--t-s)">${esc(s.biomass_id)}</span>`, `${fmt.format(nBiomass)} components`)}
    </div>

    <h3>Growth rate anchors</h3>
    <p class="sub">Anchor interval per medium axis as recorded in the model release, in h<sup>-1</sup>.</p>
    <div class="tablewrap" style="max-width:480px">
      <table class="data growth-table">
        <thead><tr><th scope="col">Medium axis</th><th scope="col">mu min</th><th scope="col">mu max</th></tr></thead>
        <tbody>${growthRows || '<tr><td colspan="3">No growth anchors recorded for this model.</td></tr>'}</tbody>
      </table>
    </div>

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
  const CAP = 300;
  const slice = list.slice(0, CAP);
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
