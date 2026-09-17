// Metabolic Atlas: main controller. One connected workflow over a persistent
// context (substrate -> product, GEM, medium): 1 Discover (pathway search +
// union map), 2 Model (GEM dashboard + browser), 3 Simulate (flux feasibility),
// 4 Engineer (strain-design analyses). Absent values render as absent; every
// count carries its denominator.

import { loadIndex, loadGraph, loadGraphMeta, loadMetIndex, fmt, downloadBlob, csvEscape } from './data.js';
import { searchPathways, carriersBySpecies, popcount } from './search.js';
import { initGemView, openReactionInBrowser, setIndexData } from './gem.js';
import { getContext, setContext, onContext, initContextBar } from './context.js';
import { GROUP_COLORS, chartBlock, hBars, histogram, heatmap, presenceMatrix } from './charts.js';

const esc = (s) => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const $ = (sel) => document.querySelector(sel);

// ---- species accent pairs (accent = graphics, ink = text-safe, wash = tint)
const SPECIES_TOKENS = {
  'Parageobacillus thermoglucosidasius': { accent: '#C0793A', ink: '#8A5222', wash: '#F6EDE1' },
  'Pseudomonas putida': { accent: '#2E6E8E', ink: '#245C77', wash: '#E8F0F4' },
  'Cupriavidus necator': { accent: '#3E8E6E', ink: '#2C6B50', wash: '#E8F2ED' },
  'Eubacterium limosum': { accent: '#7A5EA6', ink: '#5F4390', wash: '#EFEAF6' },
};
const DEFAULT_TOKENS = { accent: '#2E6E8E', ink: '#245C77', wash: '#EAF1F5' };

function setAccent(speciesName) {
  const t = SPECIES_TOKENS[speciesName] || DEFAULT_TOKENS;
  const r = document.documentElement.style;
  r.setProperty('--accent', t.accent);
  r.setProperty('--accent-ink', t.ink);
  r.setProperty('--accent-wash', t.wash);
}

// ---- app state
let INDEX = null, GRAPH = null, META = null;
let MET_INDEX = null;                    // {mid: {name, kegg, chebi}} or null when unavailable
let metList = [];                        // [{mid, name, kegg, comp, group, cur}]
let map = null;                          // map3d api or null
let selectedSpecies = new Set();
let lastResults = null;
let MODE2 = null;                        // lazily imported mode2 module api or null
let mode2Loading = false;
let highlightedIdx = null;               // pathway index currently on the map

// Display name for a metabolite id: standard name when the index carries one,
// otherwise the id itself.
function metName(mid) {
  const e = MET_INDEX && MET_INDEX[mid];
  if (e && e.name && e.name !== mid) return e.name;
  const g = GRAPH && GRAPH.metabolites[mid];
  if (g && g.n && g.n !== mid) return g.n;
  return '';
}
function metLabelHTML(mid) {
  const nm = metName(mid);
  return nm ? `${esc(nm)} <span class="mono">${esc(mid)}</span>` : `<span class="mono">${esc(mid)}</span>`;
}

// ================= boot =================
(async function boot() {
  route();
  window.addEventListener('hashchange', route);
  $('#sim-run').addEventListener('click', () => {
    if (MODE2 && lastResults) MODE2.runFeasibility(lastResults);
  });
  // endpoints set elsewhere (the Engineer pickers) reflect back into the
  // Discover inputs, so the two stages never show different context values
  onContext((c, changed) => {
    for (const kind of ['sub', 'prod']) {
      if (changed.includes(kind) && c[kind] !== pickerState[kind].mid) {
        pickerState[kind].mid = c[kind];
        const input = $(`#${kind}-input`);
        if (input) input.value = c[kind] ? (metName(c[kind]) ? `${metName(c[kind])} (${c[kind]})` : c[kind]) : '';
        if (MODE2) MODE2.onEndpointsChange(pickerState.sub.mid, pickerState.prod.mid);
      }
    }
    refreshSimUI();
  });

  try {
    INDEX = await loadIndex();
  } catch (e) {
    $('#dataset-line').textContent = `Could not load the GEM index (${e.message}). Reload the page to retry.`;
    return;
  }
  setIndexData(INDEX);
  initGemView($('#view-model'), INDEX, { onAccentChange: setAccent });
  buildSpeciesChips();
  initContextBar($('#context-bar'), {
    metLabel: (mid) => metName(mid) || mid,
    gemLabel: (acc) => acc,
  });

  try {
    [GRAPH, META] = await Promise.all([loadGraph(), loadGraphMeta()]);
  } catch (e) {
    $('#dataset-line').textContent = `Could not load the union map (${e.message}). Reload the page to retry; the Model stage works without it.`;
    $('#map-status').textContent = 'Union map unavailable. Reload the page to retry.';
    return;
  }

  try {
    MET_INDEX = await loadMetIndex();
  } catch (e) {
    MET_INDEX = null;   // pickers fall back to id-only matching and say so
  }

  metList = Object.entries(GRAPH.metabolites).map(([mid, m]) => {
    const x = MET_INDEX && MET_INDEX[mid];
    const name = (x && x.name && x.name !== mid) ? x.name : (m.n !== mid ? m.n : '');
    return { mid, name, kegg: (x && x.kegg) || '', comp: m.c, group: m.g, cur: !!m.cur };
  });
  const nCur = metList.filter(m => m.cur).length;
  $('#dataset-line').textContent =
    `${INDEX.gems.length} GEMs · ${INDEX.species.length} species · union map: ` +
    `${fmt.format(GRAPH.n_metabolites)} metabolites (${fmt.format(nCur)} currency), ` +
    `${fmt.format(GRAPH.n_reactions)} reactions, ${GRAPH.groups.length} pathway groups.`;
  $('#footer-counts').textContent =
    `Union map: ${fmt.format(GRAPH.n_metabolites)} metabolites, ${fmt.format(GRAPH.n_reactions)} reactions across ${INDEX.gems.length} GEMs.`;

  setupPicker('sub');
  setupPicker('prod');
  $('#sub-input').disabled = false;
  $('#prod-input').disabled = false;
  $('#run-search').disabled = false;
  buildExamples();
  $('#search-form').addEventListener('submit', (e) => { e.preventDefault(); runSearch(); });
  refreshSimUI();
  renderAtlasOverview();

  initMap();
})();

// ================= routing: 4 stages + legacy aliases =================
function currentStage() {
  const stageOfHash = [
    ['#/discover', 'discover'], ['#/search', 'discover'],
    ['#/model', 'model'], ['#/gems', 'model'],
    ['#/simulate', 'simulate'],
    ['#/engineer', 'engineer'], ['#/analysis', 'engineer'],
  ];
  const h = location.hash || '#/discover';
  for (const [prefix, stage] of stageOfHash) if (h.startsWith(prefix)) return stage;
  return 'discover';
}

function route() {
  const stage = currentStage();
  for (const s of ['discover', 'model', 'simulate', 'engineer']) {
    $(`#view-${s}`).hidden = s !== stage;
    $(`#nav-${s}`).setAttribute('aria-current', s === stage ? 'page' : 'false');
  }
  dockWorkbench(stage);
  if (stage === 'engineer') activateAnalysis();
  if (stage === 'simulate') activateSimulate();
}

// The map pane and the pathway-results panel are one pair of live DOM nodes,
// docked into Discover or Simulate depending on the stage, so highlights,
// feasibility chips and listeners survive the stage change.
function dockWorkbench(stage) {
  const wb = $('#workbench');
  if (!wb) return;
  if (stage === 'simulate') {
    $('#sim-dock').appendChild(wb);
  } else if (stage === 'discover' && wb.parentElement !== $('#view-discover')) {
    $('#view-discover').insertBefore(wb, $('#atlas-overview'));
  }
}

// ================= Engineer view (lazy) =================
let ANALYSIS = null, analysisLoading = false;
async function activateAnalysis() {
  if (ANALYSIS) { ANALYSIS.setActive(true); return; }
  if (analysisLoading) return;
  analysisLoading = true;
  const section = $('#view-engineer');
  try {
    for (let i = 0; i < 40 && !INDEX; i++) await new Promise(r => setTimeout(r, 200));   // boot may still be fetching the index
    if (!INDEX) {
      section.innerHTML = '<h1>Engineer</h1><p class="status error">The GEM index has not loaded; reload the page and open Engineer again.</p>';
      return;
    }
    const mod = await import('./analysis.js');
    ANALYSIS = await mod.initAnalysis(section, {
      index: INDEX,
      metName,
      setAccent,
      getEndpoints: () => ({ sub: getContext().sub, prod: getContext().prod }),
      mapAvailable: () => !!map,
      showReactionsOnMap: (rids) => {
        if (!map) return null;
        const res = map.highlightReactions(rids, accentInk());
        location.hash = '#/discover';
        $('#clear-highlight').hidden = false;
        return res;
      },
    });
    ANALYSIS.setActive(true);
  } catch (e) {
    section.innerHTML = `<h1>Engineer</h1>
      <p class="status error">The analysis module failed to load (${esc(e.message)}). Reload the page to retry; Discover and Model keep working.</p>`;
  } finally {
    analysisLoading = false;
  }
}

// ================= species chips =================
function buildSpeciesChips() {
  const row = $('#species-chips');
  row.innerHTML = '';
  for (const sp of INDEX.species) {
    selectedSpecies.add(sp.name);
    const b = document.createElement('button');
    b.type = 'button';
    b.className = 'spchip';
    b.style.setProperty('--sw', sp.accent);
    b.setAttribute('aria-pressed', 'true');
    b.dataset.name = sp.name;
    b.innerHTML = `<span class="swatch" aria-hidden="true"></span>${esc(sp.name.split(' ')[0][0])}. ${esc(sp.name.split(' ').slice(1).join(' '))} (${sp.n})`;
    b.addEventListener('click', () => {
      const on = b.getAttribute('aria-pressed') === 'true';
      if (on && selectedSpecies.size === 1) return;   // keep at least one species selected
      b.setAttribute('aria-pressed', String(!on));
      if (on) selectedSpecies.delete(sp.name); else selectedSpecies.add(sp.name);
      setAccent(selectedSpecies.size === 1 ? [...selectedSpecies][0] : null);
    });
    row.appendChild(b);
  }
}

// ================= metabolite pickers =================
const pickerState = { sub: { mid: null }, prod: { mid: null } };

function setupPicker(kind) {
  const input = $(`#${kind}-input`);
  const listbox = $(`#${kind}-listbox`);
  let options = [];        // current rendered options
  let active = -1;

  const close = () => { listbox.hidden = true; input.setAttribute('aria-expanded', 'false'); active = -1; };
  const open = () => { listbox.hidden = false; input.setAttribute('aria-expanded', 'true'); };

  // Match on id, standard name and KEGG id. Rank: exact match (0), prefix (1),
  // substring (2); ties break on shorter label.
  function scoreMatch(m, query) {
    const idL = m.mid.toLowerCase(), nmL = m.name.toLowerCase(), kgL = m.kegg.toLowerCase();
    if (idL === query || nmL === query || kgL === query) return 0;
    if (idL.startsWith(query) || nmL.startsWith(query) || kgL.startsWith(query)) return 1;
    if (idL.includes(query) || nmL.includes(query) || kgL.includes(query)) return 2;
    return -1;
  }

  function renderOptions(q) {
    const query = q.trim().toLowerCase();
    if (!query) { close(); return; }
    const scored = [];
    for (const m of metList) {
      const s = scoreMatch(m, query);
      if (s >= 0) scored.push([s, m]);
    }
    scored.sort((a, b) => a[0] - b[0]
      || (a[1].name || a[1].mid).length - (b[1].name || b[1].mid).length
      || a[1].mid.localeCompare(b[1].mid));
    const all = scored.map(x => x[1]);
    options = all.slice(0, 50);
    if (!all.length) {
      const scope = MET_INDEX
        ? `${fmt.format(metList.length)} metabolites (id, name or KEGG id)`
        : `${fmt.format(metList.length)} metabolite ids (the name index did not load; reload to retry name and KEGG matching)`;
      listbox.innerHTML = `<li class="mcap" role="presentation">No metabolite matches "${esc(q)}" among ${scope}.</li>`;
      open(); return;
    }
    const cap = all.length > options.length
      ? `<li class="mcap" role="presentation">Showing ${options.length} of ${fmt.format(all.length)} matches; keep typing to narrow.</li>` : '';
    listbox.innerHTML = options.map((m, i) => `
      <li id="${kind}-opt-${i}" role="option" aria-selected="false" data-mid="${esc(m.mid)}">
        <span class="mname-primary">${esc(m.name || m.mid)}</span>
        ${m.name ? `<span class="mid">${esc(m.mid)}</span>` : ''}
        ${m.kegg ? `<span class="mkegg">${esc(m.kegg)}</span>` : ''}
        ${m.cur ? '<span class="mname">currency</span>' : ''}
        <span class="comp-badge">${esc(m.comp)}</span>
      </li>`).join('') + cap;
    listbox.querySelectorAll('li[role="option"]').forEach(li => {
      li.addEventListener('mousedown', (e) => { e.preventDefault(); choose(li.dataset.mid); });
    });
    open();
  }

  function choose(mid) {
    pickerState[kind].mid = mid;
    input.value = metName(mid) ? `${metName(mid)} (${mid})` : mid;
    close();
    setContext(kind === 'sub' ? { sub: mid } : { prod: mid });
    if (MODE2) MODE2.onEndpointsChange(pickerState.sub.mid, pickerState.prod.mid);
  }

  function setActive(i) {
    const lis = listbox.querySelectorAll('li[role="option"]');
    if (!lis.length) return;
    active = (i + lis.length) % lis.length;
    lis.forEach((li, j) => li.setAttribute('aria-selected', String(j === active)));
    input.setAttribute('aria-activedescendant', `${kind}-opt-${active}`);
    lis[active].scrollIntoView({ block: 'nearest' });
  }

  input.addEventListener('input', () => { pickerState[kind].mid = null; renderOptions(input.value); });
  input.addEventListener('focus', () => { if (input.value && !pickerState[kind].mid) renderOptions(input.value); });
  input.addEventListener('blur', () => setTimeout(close, 120));
  input.addEventListener('keydown', (e) => {
    if (listbox.hidden && (e.key === 'ArrowDown' || e.key === 'ArrowUp')) { renderOptions(input.value); return; }
    if (e.key === 'ArrowDown') { e.preventDefault(); setActive(active + 1); }
    else if (e.key === 'ArrowUp') { e.preventDefault(); setActive(active - 1); }
    else if (e.key === 'Enter') {
      if (!listbox.hidden && active >= 0 && options[active]) { e.preventDefault(); choose(options[active].mid); }
    } else if (e.key === 'Escape') close();
  });
}

function setPickerValue(kind, mid) {
  pickerState[kind].mid = mid;
  $(`#${kind}-input`).value = metName(mid) ? `${metName(mid)} (${mid})` : mid;
  setContext(kind === 'sub' ? { sub: mid } : { prod: mid });
  if (MODE2) MODE2.onEndpointsChange(pickerState.sub.mid, pickerState.prod.mid);
}

function buildExamples() {
  const examples = [
    ['glc__D_e', 'etoh_e', 'glucose to ethanol'],
    ['glc__D_e', 'succ_e', 'glucose to succinate'],
    ['ac_e', 'etoh_e', 'acetate to ethanol'],
  ].filter(([a, b]) => GRAPH.metabolites[a] && GRAPH.metabolites[b]);
  const row = $('#example-row');
  row.innerHTML = '';
  for (const [a, b, label] of examples) {
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'btn small';
    btn.textContent = `Try: ${label}`;
    btn.addEventListener('click', () => {
      setPickerValue('sub', a);
      setPickerValue('prod', b);
      runSearch();
    });
    row.appendChild(btn);
  }
}

// ================= search (Mode 1 enumeration; feasibility is stage 3) =================
function runSearch() {
  const sub = pickerState.sub.mid, prod = pickerState.prod.mid;
  const summary = $('#results-summary');
  if (!GRAPH || !META) return;
  if (!sub || !prod) {
    summary.innerHTML = `<span class="status">Pick both a substrate and a product from the list (type a name, BiGG id or KEGG id to search ${fmt.format(metList.length)} metabolites).</span>`;
    return;
  }
  if (sub === prod) {
    summary.innerHTML = '<span class="status">Substrate and product are the same metabolite; nothing to search.</span>';
    return;
  }
  summary.textContent = 'Searching…';
  setTimeout(() => {
    const t0 = performance.now();
    const res = searchPathways(GRAPH, META, sub, prod, [...selectedSpecies]);
    const ms = Math.round(performance.now() - t0);
    lastResults = { res, sub, prod, species: [...selectedSpecies] };
    renderResults(res, sub, prod, ms);
    if (res.pathways.length && map) showPathwayOnMap(res.pathways[0], 0, { animate: true });
    else clearMapHighlight();
    refreshSimUI();
  }, 30);
}

function terminationLine(t, nFound) {
  if (t.complete) {
    return `<p class="termination">Search complete: all simple pathways up to ${t.maxDepth} steps enumerated ` +
      `(${fmt.format(t.expansions)} edge expansions). Longer pathways are not enumerated.</p>`;
  }
  const parts = [];
  if (t.hitPaths) parts.push(`stopped at the ${t.maxPaths}-pathway cap`);
  if (t.hitBudget) parts.push(`stopped at the ${fmt.format(t.maxExpansions)}-expansion budget`);
  return `<p class="termination capped">Search incomplete: ${parts.join('; ')} after ${fmt.format(t.expansions)} edge expansions, ` +
    `within the ${t.maxDepth}-step limit. ${nFound ? 'Further pathways may exist beyond these caps.' : 'A pathway may exist beyond these caps.'}</p>`;
}

function renderResults(res, sub, prod, ms) {
  const { pathways, termination } = res;
  const summary = $('#results-summary');
  const body = $('#results-body');
  const nGems = META.accs.length;
  const spNote = lastResults.species.length < INDEX.species.length
    ? ` · species filter: ${lastResults.species.length} of ${INDEX.species.length} species` : '';

  if (!pathways.length) {
    summary.innerHTML = `<div class="result-summary"><strong>No pathway found</strong> from
      ${metLabelHTML(sub)} to ${metLabelHTML(prod)}
      within ${termination.maxDepth} steps for the selected species${spNote} (${ms} ms).</div>`;
    body.innerHTML = `<div class="card empty-state">
      ${terminationLine(termination, 0)}
      <p>Currency metabolites are never used as intermediates, so routes that pass only through
      ATP, NAD(H), water and similar carriers are not enumerated. Compartment matters:
      <span class="mono">_e</span>, <span class="mono">_p</span> and <span class="mono">_c</span>
      forms of a metabolite are distinct nodes joined by transport reactions.</p></div>`;
    return;
  }

  const best = pathways[0];
  summary.innerHTML = `<div class="result-summary">
    <strong>${fmt.format(pathways.length)} pathway${pathways.length > 1 ? 's' : ''}</strong> from
    ${metLabelHTML(sub)} to ${metLabelHTML(prod)} ·
    shortest ${best.len} step${best.len > 1 ? 's' : ''} ·
    best carried end-to-end by ${best.carriers} of ${nGems} GEMs${spNote} · ${ms} ms</div>`;

  body.innerHTML = '';
  const CARD_BATCH = 20;
  let shown = 0;
  const showMore = () => {
    const next = pathways.slice(shown, shown + CARD_BATCH);
    next.forEach((pw, i) => body.insertBefore(pathwayCard(pw, shown + i), moreBar));
    shown += next.length;
    moreBar.querySelector('span').textContent = `Showing ${shown} of ${fmt.format(pathways.length)} pathways.`;
    moreBar.querySelector('button').hidden = shown >= pathways.length;
  };
  const moreBar = document.createElement('div');
  moreBar.className = 'card';
  moreBar.innerHTML = `<span class="status"></span> <button class="btn small" type="button">Show ${CARD_BATCH} more</button>`;
  moreBar.querySelector('button').addEventListener('click', showMore);
  body.appendChild(moreBar);
  showMore();

  const term = document.createElement('div');
  term.innerHTML = terminationLine(termination, pathways.length);
  body.appendChild(term);

  body.appendChild(searchCharts(pathways));
}

// Two designed charts over the current result set: the pathway-length
// distribution and the strain-presence overview matrix.
function searchCharts(pathways) {
  const wrap = document.createElement('div');
  wrap.className = 'card';

  const byLen = new Map();
  for (const pw of pathways) byLen.set(pw.len, (byLen.get(pw.len) || 0) + 1);
  const lens = [...byLen.keys()].sort((a, b) => a - b);
  const bins = [];
  for (let L = lens[0]; L <= lens[lens.length - 1]; L++) bins.push({ x: L, count: byLen.get(L) || 0 });
  const lenChart = chartBlock(
    `Pathway length distribution (${fmt.format(pathways.length)} pathways)`,
    histogram(bins, { width: 330, height: 120, xLabel: 'steps' }));

  const nRows = Math.min(pathways.length, 15);
  const rows = pathways.slice(0, nRows);
  const rowLabels = rows.map((pw, i) => `#${i + 1} · ${pw.len} step${pw.len > 1 ? 's' : ''}`);
  const colLabels = META.accs.map(a => `${a.acc} (${a.sp})`);
  const colColors = META.accs.map(a => (SPECIES_TOKENS[a.sp] || DEFAULT_TOKENS).accent);
  const cells = rows.map(pw => META.accs.map((a, j) => ((pw.mask >> BigInt(j)) & 1n) === 1n));
  const legend = INDEX.species.map(sp =>
    `<span class="ch-lg"><span class="ch-sw" style="background:${sp.accent}"></span>${esc(sp.name)} (${sp.n})</span>`).join('');
  const presChart = chartBlock(
    `Strain presence per pathway (${nRows} of ${fmt.format(pathways.length)} pathways shown × ${META.accs.length} GEMs)`,
    `<div class="chartwrap">${presenceMatrix(rowLabels, colLabels, cells, colColors, { cellW: 12, labelW: 86 })}</div>
     <div class="ch-legend">${legend}</div>`,
    'A filled cell: the GEM contains at least one reaction for every step. Flux feasibility under a medium is the Simulate stage.');

  wrap.innerHTML = lenChart + presChart;
  return wrap;
}

function pathwayCard(pw, idx) {
  const nGems = META.accs.length;
  const d = document.createElement('details');
  d.className = 'pcard';
  d.dataset.pwIdx = idx;
  const chain = pw.mets.map(m => metName(m) || m).join(' → ');
  d.innerHTML = `
    <summary>
      <span class="rank">#${idx + 1}</span>
      <span class="plen">${pw.len} step${pw.len > 1 ? 's' : ''}</span>
      <span class="carriers">${pw.carriers} of ${nGems} GEMs carry every step</span>
      <span class="feas-slot"></span>
      <span class="chain" title="${esc(pw.mets.join(' → '))}">${esc(chain)}</span>
    </summary>
    <div class="pcard-body"></div>`;
  const bodyEl = d.querySelector('.pcard-body');

  // Live map highlight: hover or keyboard focus previews, opening the card
  // animates the walk. The map is optional; without it these are no-ops.
  const preview = () => { if (map && highlightedIdx !== idx) showPathwayOnMap(pw, idx, { animate: false }); };
  if (window.matchMedia('(hover: hover) and (pointer: fine)').matches) {
    d.querySelector('summary').addEventListener('mouseenter', preview);
  }
  d.querySelector('summary').addEventListener('focus', preview);
  d.addEventListener('toggle', () => { if (d.open && map) showPathwayOnMap(pw, idx, { animate: true }); });

  const stepsEl = document.createElement('div');
  pw.steps.forEach((st, i) => {
    const altMax = 6;
    const shown = st.rxns.slice(0, altMax);
    const stepDiv = document.createElement('div');
    stepDiv.className = 'step';
    stepDiv.dataset.step = i;
    const fromNm = metName(st.from), toNm = metName(st.to);
    stepDiv.innerHTML = `<div class="step-mets">Step ${i + 1}: ${fromNm ? esc(fromNm) + ' ' : ''}<span class="mono">${esc(st.from)}</span> → ${toNm ? esc(toNm) + ' ' : ''}<span class="mono">${esc(st.to)}</span></div>` +
      shown.map(r => `
        <div class="rxn-alt" data-rxn="${esc(r.id)}">
          <button class="rid btn small" type="button" data-rid="${esc(r.id)}" title="Open in the Model stage">${esc(r.id)}</button>
          ${r.dir === 'rev' ? '<span class="ec">(reverse of written direction)</span>' : ''}
          ${r.ec.length ? `<span class="ec">EC ${esc(r.ec.join(', '))}</span>` : ''}
          <span class="rn">${esc(r.name || '')}</span>
          <span class="flux-slot"></span>
        </div>`).join('') +
      (st.rxns.length > altMax ? `<p class="termination">Showing ${altMax} of ${st.rxns.length} alternative reactions for this step.</p>` : '');
    stepsEl.appendChild(stepDiv);
  });
  bodyEl.appendChild(stepsEl);

  // presence matrix
  const matrix = document.createElement('div');
  matrix.className = 'matrix';
  const bySp = carriersBySpecies(pw.mask, META);
  for (const [sp, rec] of bySp) {
    const t = SPECIES_TOKENS[sp] || DEFAULT_TOKENS;
    const row = document.createElement('div');
    row.className = 'sprow';
    row.innerHTML = `<span class="spname">${esc(sp)}</span>
      <span class="spcount">${rec.present} of ${rec.total} strains</span>
      <span class="cells">${rec.cells.map(c =>
        `<span class="cell ${c.on ? 'on' : 'off'}" ${c.on ? `style="background:${t.accent}"` : ''}
          title="${esc(c.acc)}: ${c.on ? 'carries every step' : 'missing at least one step'}"></span>`).join('')}</span>`;
    matrix.appendChild(row);
  }
  const note = document.createElement('p');
  note.className = 'termination';
  note.textContent = 'A filled cell means the strain’s GEM contains at least one reaction for every step. Whether the strain’s own bounds and a medium permit flux through every step is a Simulate check.';
  matrix.appendChild(note);
  bodyEl.appendChild(matrix);

  // Feasibility detail slot (filled by the Simulate stage when a run covers this card)
  const m2slot = document.createElement('div');
  m2slot.className = 'mode2-slot';
  bodyEl.appendChild(m2slot);

  // actions
  const actions = document.createElement('div');
  actions.className = 'cardactions';
  const btnMap = document.createElement('button');
  btnMap.className = 'btn small'; btnMap.type = 'button';
  btnMap.textContent = map ? 'Show on 3D map' : '3D map unavailable';
  btnMap.disabled = !map;
  btnMap.addEventListener('click', () => {
    if (!map) return;
    showPathwayOnMap(pw, idx, { animate: true });
    $('#map-pane').scrollIntoView({ behavior: 'smooth', block: 'nearest' });
  });
  const btnCsv = document.createElement('button');
  btnCsv.className = 'btn small'; btnCsv.type = 'button'; btnCsv.textContent = 'Export CSV';
  btnCsv.addEventListener('click', () => exportPathway(pw, idx, 'csv'));
  const btnJson = document.createElement('button');
  btnJson.className = 'btn small'; btnJson.type = 'button'; btnJson.textContent = 'Export JSON';
  btnJson.addEventListener('click', () => exportPathway(pw, idx, 'json'));
  actions.append(btnMap, btnCsv, btnJson);
  bodyEl.appendChild(actions);

  bodyEl.querySelectorAll('button.rid').forEach(b => b.addEventListener('click', () => {
    openReactionInBrowser(b.dataset.rid);
    location.hash = '#/model';
  }));
  return d;
}

function exportPathway(pw, idx, kind) {
  const presence = META.accs.map((a, i) => ({
    acc: a.acc, species: a.sp, carries_every_step: ((pw.mask >> BigInt(i)) & 1n) === 1n,
  }));
  const meta = {
    substrate: lastResults.sub, product: lastResults.prod,
    species_filter: lastResults.species, rank: idx + 1, steps: pw.len,
    carriers: `${pw.carriers} of ${META.accs.length} GEMs`,
    note: 'Enzyme-presence pathway; flux feasibility under a medium is a Simulate-stage check and is not included in this export.',
  };
  if (kind === 'json') {
    downloadBlob(JSON.stringify({
      ...meta,
      metabolites: pw.mets,
      reactions_per_step: pw.steps.map(s => ({ from: s.from, to: s.to, reactions: s.rxns.map(r => ({ id: r.id, name: r.name, direction: r.dir, ec: r.ec })) })),
      presence,
    }, null, 1), `pathway_${idx + 1}_${lastResults.sub}_to_${lastResults.prod}.json`, 'application/json');
  } else {
    const lines = [
      `# pathway ${idx + 1}: ${lastResults.sub} to ${lastResults.prod}; ${pw.len} steps; carried by ${pw.carriers} of ${META.accs.length} GEMs`,
      'step,from,to,reaction_ids,directions,ec_numbers',
      ...pw.steps.map((s, i) => [i + 1, s.from, s.to, s.rxns.map(r => r.id).join(';'), s.rxns.map(r => r.dir).join(';'), s.rxns.flatMap(r => r.ec).join(';')].map(csvEscape).join(',')),
      '',
      'acc,species,carries_every_step',
      ...presence.map(p => [p.acc, p.species, p.carries_every_step].map(csvEscape).join(',')),
    ];
    downloadBlob(lines.join('\n'), `pathway_${idx + 1}_${lastResults.sub}_to_${lastResults.prod}.csv`, 'text/csv');
  }
}

// ================= Atlas overview (Discover, below the workbench) =================
function renderAtlasOverview() {
  const host = $('#atlas-overview');
  if (!host || !GRAPH || !META) return;

  // model size per species: mean over each species' GEMs, from the registry
  const bySp = new Map();
  for (const g of INDEX.gems) {
    if (!bySp.has(g.species)) bySp.set(g.species, []);
    bySp.get(g.species).push(g);
  }
  const sizeChart = (key, label) => {
    const rows = [...bySp.entries()].map(([sp, gems]) => {
      const vals = gems.map(g => g[key]);
      const mean = Math.round(vals.reduce((a, b) => a + b, 0) / vals.length);
      return {
        label: sp.split(' ')[0][0] + '. ' + sp.split(' ').slice(1).join(' '),
        value: mean,
        color: (SPECIES_TOKENS[sp] || DEFAULT_TOKENS).accent,
        note: `(n=${gems.length}, ${fmt.format(Math.min(...vals))} to ${fmt.format(Math.max(...vals))})`,
      };
    });
    return chartBlock(`${label} per GEM, mean per species (${INDEX.gems.length} GEMs)`,
      hBars(rows, { showPct: false }));
  };

  // reaction sharing: how many of the 35 GEMs carry each union reaction
  const carrierCounts = new Array(META.accs.length + 1).fill(0);
  let nWithMask = 0;
  for (const id of Object.keys(META.rxns)) {
    const m = META.rxns[id].m;
    if (m == null) continue;
    nWithMask++;
    carrierCounts[popcount(BigInt('0x' + m))]++;
  }
  // carrierCounts is dense (pre-filled with 0 for every carrier count), so an
  // index read is always a measured count, never an absent lookup.
  const nCore = carrierCounts.slice(34).reduce((a, b) => a + b, 0);
  const nSingle = carrierCounts[1];
  const shareBins = [];
  for (let k = 1; k <= META.accs.length; k++) {
    shareBins.push({
      x: k, count: carrierCounts[k],
      color: k >= 34 ? 'var(--accent)' : (k === 1 ? '#98948C' : 'rgba(42,98,176,0.55)'),
      tick: (k === 1 || k % 5 === 0 || k === META.accs.length),
    });
  }
  const nZero = carrierCounts[0];
  const shareChart = chartBlock(
    `Reaction sharing across the ${META.accs.length} GEMs (${fmt.format(nWithMask)} union reactions)`,
    `<div class="chartwrap">${histogram(shareBins, { width: 560, height: 150, xLabel: 'number of GEMs carrying the reaction', barLabel: (b) => b.x === 1 || b.x >= 34 })}</div>`,
    `Core (in 34 or 35 GEMs): ${fmt.format(nCore)} of ${fmt.format(nWithMask)} · in a single GEM: ${fmt.format(nSingle)} of ${fmt.format(nWithMask)}.` +
    (nZero ? ` ${fmt.format(nZero)} of ${fmt.format(nWithMask)} union reactions are carried by no individual GEM in this release and are not drawn (methods page, union map).` : ''));

  // pathway-group coverage: species x group, union reactions carried by >= 1 GEM
  const groups = GRAPH.groups;
  const spNames = INDEX.species.map(s => s.name);
  const counts = spNames.map(() => groups.map(() => 0));
  const groupTotals = groups.map(() => 0);
  for (const r of GRAPH.reactions) {
    const gi = groups.indexOf(r.g);
    if (gi < 0) continue;
    groupTotals[gi]++;
    for (const sp of r.sp) {
      const si = spNames.indexOf(sp);
      if (si >= 0) counts[si][gi]++;
    }
  }
  const heatChart = chartBlock(
    `Pathway-group coverage by species (union reactions per group, ${fmt.format(GRAPH.n_reactions)} total)`,
    `<div class="chartwrap">${heatmap(
      spNames.map(sp => sp.split(' ')[0][0] + '. ' + sp.split(' ').slice(1).join(' ')),
      groups.map((g, i) => `${g} (${fmt.format(groupTotals[i])})`),
      counts, { cellW: 52, labelW: 150 })}</div>`,
    'A cell counts the union reactions of the group carried by at least one GEM of the species; the column header carries the group’s union total.');

  host.innerHTML = `
    <h2>Atlas overview</h2>
    <p class="sub">The 35 models and the union network they share, before any search.</p>
    <div class="chart-grid">
      ${sizeChart('reactions', 'Reactions')}
      ${sizeChart('genes', 'Genes')}
      ${sizeChart('metabolites', 'Metabolites')}
    </div>
    <div class="chart-grid">
      ${shareChart}
      ${heatChart}
    </div>`;
}

// ================= map highlight orchestration =================
function accentInk() {
  return getComputedStyle(document.documentElement).getPropertyValue('--accent-ink').trim() || '#245C77';
}

// weights: optional array (one per step) of |flux| magnitudes; the map scales
// edge thickness by them. animate steps the walk substrate to product.
function showPathwayOnMap(pw, idx, { animate = false, weights = null } = {}) {
  if (!map) return;
  highlightedIdx = idx;
  map.highlightPathway(pw.mets, accentInk(), { animate, weights });
  $('#clear-highlight').hidden = false;
}

function clearMapHighlight() {
  if (map) { map.clearHighlight(); map.clearReactionHighlight(); }
  highlightedIdx = null;
  $('#clear-highlight').hidden = true;
}

// ================= Simulate stage (lazy flux engine) =================
async function activateSimulate() {
  refreshSimUI();
  if (MODE2 || mode2Loading) return;
  mode2Loading = true;
  const hostEl = $('#mode2-host');
  try {
    const mod = await import('./mode2.js');
    MODE2 = await mod.initMode2(hostEl, {
      index: INDEX, graphMeta: META,
      metName, metLabelHTML, setAccent,
      getEndpoints: () => ({ sub: pickerState.sub.mid, prod: pickerState.prod.mid }),
      findCard: (i) => document.querySelector(`.pcard[data-pw-idx="${i}"]`),
      showFluxOnMap: (pw, i, weights) => showPathwayOnMap(pw, i, { animate: false, weights }),
      resultsSummary: () => $('#results-summary'),
      speciesTokens: SPECIES_TOKENS,
      onStateChange: () => refreshSimUI(),
    });
    MODE2.onEndpointsChange(pickerState.sub.mid, pickerState.prod.mid);
  } catch (e) {
    hostEl.innerHTML = `<div class="card"><p class="status error">The flux engine failed to load (${esc(e.message)}).
      Reload the page to retry. Discover and Model keep working.</p></div>`;
  } finally {
    mode2Loading = false;
    refreshSimUI();
  }
}

function refreshSimUI() {
  const btn = $('#sim-run'), status = $('#sim-status');
  if (!btn) return;
  const disable = (msg) => { btn.disabled = true; status.innerHTML = msg; };
  if (!GRAPH || !META) return disable('The dataset is still loading.');
  if (!lastResults || !lastResults.res.pathways.length) {
    return disable('No pathways to test yet: run a substrate to product search in <a href="#/discover">Discover</a> first.');
  }
  if (mode2Loading || !MODE2) return disable('The flux engine is loading…');
  const why = MODE2.notReadyReason();
  if (why) return disable(esc(why));
  const n = lastResults.res.pathways.length;
  btn.disabled = false;
  status.innerHTML = `${fmt.format(n)} pathway${n > 1 ? 's' : ''} from ${metLabelHTML(lastResults.sub)}
    to ${metLabelHTML(lastResults.prod)} ready to test; the shortest ${Math.min(10, n)} of ${fmt.format(n)} are tested automatically, the rest per card.`;
}

// ================= 3D map =================
async function initMap() {
  const pane = $('#map-pane');
  const statusEl = $('#map-status');
  try {
    const { createMap } = await import('./map3d.js');
    const tip = $('#map-tip');
    map = await createMap(pane, GRAPH, GROUP_COLORS, {
      onHover(hit) {
        if (!hit) { tip.hidden = true; return; }
        const nm = metName(hit.mid);
        tip.innerHTML = `<span class="tid">${esc(hit.mid)}</span><br>
          ${nm ? esc(nm) + '<br>' : ''}
          ${esc(hit.met.g)} · compartment ${esc(hit.met.c)}${hit.met.cur ? ' · currency' : ''}`;
        tip.style.left = Math.min(hit.x + 12, pane.clientWidth - 240) + 'px';
        tip.style.top = (hit.y + 12) + 'px';
        tip.hidden = false;
      },
    });
    statusEl.hidden = true;
    $('#map-hud').hidden = false;
    buildLegend();
    // currency metabolites start hidden so the three shells read as structure;
    // the HUD checkbox brings them back.
    map.setCurrencyVisible(false);
    $('#toggle-currency').checked = false;
    $('#toggle-currency').addEventListener('change', (e) => map.setCurrencyVisible(e.target.checked));
    $('#toggle-labels').addEventListener('change', (e) => map.setLabelsVisible(e.target.checked));
    $('#reset-view').addEventListener('click', () => map.resetView());
    $('#clear-highlight').addEventListener('click', clearMapHighlight);
  } catch (e) {
    map = null;
    const reason = e.message === 'webgl-unavailable'
      ? 'This browser session does not provide WebGL, which the 3D map needs.'
      : e.message === 'cdn-unavailable'
        ? 'The three.js library could not be loaded from the CDN (network unavailable?).'
        : `The 3D map failed to start (${e.message}).`;
    statusEl.innerHTML = `${esc(reason)}<br>
      The union map holds ${fmt.format(GRAPH.n_metabolites)} metabolites and
      ${fmt.format(GRAPH.n_reactions)} reactions on three compartment shells.
      Pathway search and the Model stage work without the 3D view.`;
  }
}

function buildLegend() {
  const el = $('#map-legend');
  el.hidden = false;
  el.innerHTML =
    '<span class="lg" style="font-weight:600">Shells:</span>' +
    '<span class="lg">outer = extracellular + exchange</span>' +
    '<span class="lg">middle = periplasm</span>' +
    '<span class="lg">inner = cytosol</span>' +
    '<span class="lg" style="font-weight:600;margin-left:6px">Groups:</span>' + GRAPH.groups.map(g =>
      `<span class="lg"><span class="swatch" style="background:${GROUP_COLORS[g] || '#98948C'}"></span>${esc(g)}</span>`).join('');
}
