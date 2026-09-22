// Metabolic Atlas: main controller. One connected workflow over a persistent
// context (substrate -> product, GEM, medium): 1 Discover (pathway search +
// union map), 2 Model (GEM dashboard + browser), 3 Simulate (flux feasibility),
// 4 Engineer (strain-design analyses). Absent values render as absent; every
// count carries its denominator.

import { loadIndex, loadGraph, loadGraphMeta, loadMetIndex, loadGem, fmt, downloadBlob, csvEscape, geneLabel, geneLabelHTML, subsystemLabel } from './data.js';
import { searchPathways, carriersBySpecies, popcount } from './search.js';
import { initGemView, openReactionInBrowser, setIndexData } from './gem.js';
import { getContext, setContext, onContext, initContextBar } from './context.js';
import { GROUP_COLORS, SEARCH_COLORS, chartBlock, hBars, histogram, heatmap, presenceMatrix } from './charts.js';

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
let ctxGem = null;                       // loaded GEM json of the context accession, for gene/equation labels
let fluxLayerOn = false;                 // pale-orange flux edges currently drawn

// Results-list view state (P1 controls); reset on every new search.
const viewState = { lenFilter: null, feasOnly: false, sort: 'len' };

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
  $('#sim-cancel').addEventListener('click', () => { if (MODE2) MODE2.cancelRun(); });
  // keep a loaded copy of the context GEM so pathway cards can label genes,
  // equations and subsystems; absent stays absent until the load lands
  onContext(async (c, changed) => {
    if (!changed.includes('gem')) return;
    if (!c.gem) { ctxGem = null; return; }
    try {
      const g = await loadGem(c.gem);
      if (getContext().gem === c.gem) { ctxGem = g; enrichOpenCards(); }
    } catch { ctxGem = null; }
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

  metList = Object.entries(GRAPH.metabolites)
    .filter(([, m]) => !m.core)   // biomass_core is a map visual, not a searchable metabolite
    .map(([mid, m]) => {
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
    viewState.lenFilter = null; viewState.feasOnly = false; viewState.sort = 'len';
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
    summary.innerHTML = `<div class="result-head">
      <div class="result-line"><span class="result-count none">No pathway found</span>
        <span class="result-route">${metLabelHTML(sub)} <span class="route-arrow" aria-hidden="true">→</span> ${metLabelHTML(prod)}</span></div>
      <div class="result-meta">within ${termination.maxDepth} steps for the selected species${spNote} · search ${ms} ms</div>
    </div>`;
    body.innerHTML = `<div class="card empty-state">
      ${terminationLine(termination, 0)}
      <p>Currency metabolites are never used as intermediates, so routes that pass only through
      ATP, NAD(H), water and similar carriers are not enumerated. Compartment matters:
      <span class="mono">_e</span>, <span class="mono">_p</span> and <span class="mono">_c</span>
      forms of a metabolite are distinct nodes joined by transport reactions.</p></div>`;
    return;
  }

  const best = pathways[0];
  summary.innerHTML = `<div class="result-head">
    <div class="result-line">
      <span class="result-count">${fmt.format(pathways.length)}</span>
      <span class="result-word">pathway${pathways.length > 1 ? 's' : ''}</span>
      <span class="result-route">${metLabelHTML(sub)} <span class="route-arrow" aria-hidden="true">→</span> ${metLabelHTML(prod)}</span>
    </div>
    <div class="result-meta">shortest ${best.len} step${best.len > 1 ? 's' : ''} ·
      best carried end-to-end by ${best.carriers} of ${nGems} GEMs${spNote} · search ${ms} ms</div>
  </div>`;

  body.innerHTML = '';
  const controls = document.createElement('div');
  controls.className = 'card list-controls';
  controls.id = 'list-controls';
  body.appendChild(controls);

  const listHost = document.createElement('div');
  listHost.className = 'pathway-list';
  listHost.id = 'pathway-list';
  body.appendChild(listHost);

  const term = document.createElement('div');
  term.innerHTML = terminationLine(termination, pathways.length);
  body.appendChild(term);

  body.appendChild(searchCharts(pathways));

  renderListControls();
  renderCardList();
}

// ---- results view: interactive length histogram + sort + feasible-only ----
function runView() {
  return (MODE2 && MODE2.getRun && MODE2.getRun()) || null;
}

// Rebuild the controls (histogram counts change when a feasibility run lands).
function renderListControls() {
  const host = document.getElementById('list-controls');
  if (!host || !lastResults) return;
  const pathways = lastResults.res.pathways;
  const run = runView();

  const byLen = new Map();
  for (const pw of pathways) {
    if (!byLen.has(pw.len)) byLen.set(pw.len, { total: 0, feasible: 0, infeasible: 0 });
    byLen.get(pw.len).total++;
  }
  if (run) {
    pathways.forEach((pw, i) => {
      const r = run.results.get(i);
      if (!r) return;
      if (r.testable && r.feasible) byLen.get(pw.len).feasible++;
      else byLen.get(pw.len).infeasible++;     // tested infeasible, or not carried by this GEM
    });
  }
  const lens = [...byLen.keys()].sort((a, b) => a - b);
  const maxTotal = Math.max(...lens.map(L => byLen.get(L).total), 1);
  const BAR_H = 56;
  const nFeas = run ? run.feasibleIdx.length : null;
  const nTested = run ? run.results.size : null;

  const bars = lens.map(L => {
    const b = byLen.get(L);
    const untested = b.total - b.feasible - b.infeasible;
    const h = (v) => Math.round(v / maxTotal * BAR_H);
    const pressed = viewState.lenFilter === L;
    const segs =
      (b.feasible ? `<span class="lenseg" style="height:${Math.max(h(b.feasible), 2)}px;background:${SEARCH_COLORS.best}"></span>` : '') +
      (b.infeasible ? `<span class="lenseg" style="height:${Math.max(h(b.infeasible), 2)}px;background:${SEARCH_COLORS.infeasible}"></span>` : '') +
      (untested ? `<span class="lenseg" style="height:${Math.max(h(untested), 2)}px;background:${SEARCH_COLORS.untested}"></span>` : '');
    const tip = run
      ? `${b.total} pathway${b.total > 1 ? 's' : ''} of ${L} steps: ${b.feasible} feasible, ${b.infeasible} infeasible or not carried, ${untested} untested`
      : `${b.total} pathway${b.total > 1 ? 's' : ''} of ${L} steps (feasibility not tested yet)`;
    return `<button type="button" class="lenbar" data-len="${L}" aria-pressed="${pressed}" title="${esc(tip)}" aria-label="${esc(tip)}${pressed ? '; filter active' : '; filter the list to this length'}">
      <span class="lencount">${b.total}</span>
      <span class="lenstack" aria-hidden="true">${segs}</span>
      <span class="lenlab">${L}</span>
    </button>`;
  }).join('');

  const runNote = run
    ? `feasibility on <span class="mono">${esc(run.acc)}</span> · ${esc(run.medium || '')}: ${nFeas} of ${nTested} tested feasible (${fmt.format(pathways.length)} found)`
    : `${fmt.format(pathways.length)} pathways found; feasibility not tested yet (stage 3, Simulate)`;

  host.innerHTML = `
    <div class="lc-hist">
      <div class="chart-title">Pathway lengths (${fmt.format(pathways.length)} pathways); choose a bar to filter</div>
      <div class="lenchart" role="group" aria-label="Pathway length filter">${bars}</div>
      <div class="lenlegend">
        <span class="ch-lg"><span class="ch-sw" style="background:${SEARCH_COLORS.best}"></span>feasible</span>
        <span class="ch-lg"><span class="ch-sw" style="background:${SEARCH_COLORS.infeasible}"></span>infeasible / not carried</span>
        <span class="ch-lg"><span class="ch-sw" style="background:${SEARCH_COLORS.untested}"></span>not tested</span>
        <span class="lennote">${runNote}</span>
      </div>
    </div>
    <div class="listbar">
      <div class="field">
        <label for="sort-select">Sort</label>
        <select id="sort-select">
          <option value="len" ${viewState.sort === 'len' ? 'selected' : ''}>Shortest first</option>
          <option value="flux" ${viewState.sort === 'flux' ? 'selected' : ''} ${run ? '' : 'disabled'}>Highest product flux first${run ? '' : ' (run Simulate)'}</option>
        </select>
      </div>
      <label class="mode feas-toggle" style="border-style:solid">
        <input type="checkbox" id="feas-only" ${viewState.feasOnly ? 'checked' : ''} ${run ? '' : 'disabled'}>
        Feasible only${run ? '' : ' (run Simulate first)'}
      </label>
      ${viewState.lenFilter != null ? `<button class="btn small" type="button" id="len-clear">Clear length filter (${viewState.lenFilter} steps)</button>` : ''}
    </div>`;

  host.querySelectorAll('.lenbar').forEach(b => b.addEventListener('click', () => {
    const L = +b.dataset.len;
    viewState.lenFilter = viewState.lenFilter === L ? null : L;
    renderListControls();
    renderCardList();
  }));
  host.querySelector('#sort-select').addEventListener('change', (e) => {
    viewState.sort = e.target.value;
    renderCardList();
  });
  const ft = host.querySelector('#feas-only');
  if (ft) ft.addEventListener('change', (e) => { viewState.feasOnly = e.target.checked; renderCardList(); });
  const lc = host.querySelector('#len-clear');
  if (lc) lc.addEventListener('click', () => { viewState.lenFilter = null; renderListControls(); renderCardList(); });
}

// The filtered, sorted card list; mode2's observer re-decorates every card it
// renders with the stored feasibility, FVA and sampling results.
function renderCardList() {
  const listHost = document.getElementById('pathway-list');
  if (!listHost || !lastResults) return;
  const pathways = lastResults.res.pathways;
  const run = runView();

  let view = pathways.map((pw, idx) => ({ pw, idx }));
  if (viewState.lenFilter != null) view = view.filter(v => v.pw.len === viewState.lenFilter);
  if (viewState.feasOnly && run) view = view.filter(v => {
    const r = run.results.get(v.idx);
    return r && r.testable && r.feasible;
  });
  if (viewState.sort === 'flux' && run) {
    const flux = (i) => {
      const r = run.results.get(i);
      return (r && r.testable && r.feasible && r.productFlux != null) ? r.productFlux : -Infinity;
    };
    view.sort((a, b) => flux(b.idx) - flux(a.idx) || a.pw.len - b.pw.len || a.idx - b.idx);
  }

  listHost.innerHTML = '';
  const filtered = view.length !== pathways.length;
  const moreBar = document.createElement('div');
  moreBar.className = 'card morebar';
  moreBar.innerHTML = `<span class="status"></span> <button class="btn small" type="button"></button>`;
  const CARD_BATCH = 20;
  let shown = 0;
  const showMore = () => {
    const next = view.slice(shown, shown + CARD_BATCH);
    next.forEach(v => listHost.insertBefore(pathwayCard(v.pw, v.idx), moreBar));
    shown += next.length;
    moreBar.querySelector('span').textContent =
      `Showing ${shown} of ${fmt.format(view.length)}${filtered ? ` matching pathways (of ${fmt.format(pathways.length)} found)` : ' pathways'}.`;
    const btn = moreBar.querySelector('button');
    btn.textContent = `Show ${Math.min(CARD_BATCH, view.length - shown)} more`;
    btn.hidden = shown >= view.length;
  };
  moreBar.querySelector('button').addEventListener('click', showMore);
  listHost.appendChild(moreBar);
  if (!view.length) {
    moreBar.querySelector('span').textContent = viewState.feasOnly
      ? `No feasible pathway matches the current filter (0 of ${fmt.format(pathways.length)} found pathways).`
      : `No pathway matches the current filter (0 of ${fmt.format(pathways.length)} found pathways).`;
    moreBar.querySelector('button').hidden = true;
  } else {
    showMore();
  }
}

// The strain-presence overview matrix over the current result set.
function searchCharts(pathways) {
  const wrap = document.createElement('div');
  wrap.className = 'card';

  const nRows = Math.min(pathways.length, 15);
  const rows = pathways.slice(0, nRows);
  const rowLabels = rows.map((pw, i) => `#${i + 1} · ${pw.len} step${pw.len > 1 ? 's' : ''}`);
  const colLabels = META.accs.map(a => `${a.acc} (${a.sp})`);
  const colColors = META.accs.map(a => (SPECIES_TOKENS[a.sp] || DEFAULT_TOKENS).accent);
  const cells = rows.map(pw => META.accs.map((a, j) => ((pw.mask >> BigInt(j)) & 1n) === 1n));
  const legend = INDEX.species.map(sp =>
    `<span class="ch-lg"><span class="ch-sw" style="background:${sp.accent}"></span>${esc(sp.name)} (${sp.n})</span>`).join('');
  wrap.innerHTML = chartBlock(
    `Strain presence per pathway (${nRows} of ${fmt.format(pathways.length)} pathways shown × ${META.accs.length} GEMs)`,
    `<div class="chartwrap">${presenceMatrix(rowLabels, colLabels, cells, colColors, { cellW: 12, labelW: 86 })}</div>
     <div class="ch-legend">${legend}</div>`,
    'A filled cell: the GEM contains at least one reaction for every step. Flux feasibility under a medium is the Simulate stage.');
  return wrap;
}

// Equation string for a reaction of the context GEM, bounds-aware arrow.
function gemEquation(g, rid) {
  const r = g.reactions.find(x => x.id === rid);
  if (!r) return null;
  const lhs = [], rhs = [];
  for (const [m, c] of Object.entries(r.stoich)) {
    const coef = Math.abs(c) === 1 ? '' : `${Math.round(Math.abs(c) * 1e4) / 1e4} `;
    (c < 0 ? lhs : rhs).push(`${coef}${m}`);
  }
  const arrow = (r.lb < 0 && r.ub > 0) ? '⇌' : (r.ub <= 0 && r.lb < 0) ? '←' : '→';
  return `${lhs.join(' + ')} ${arrow} ${rhs.join(' + ')}`;
}

// Typeset equation for a reaction of the context GEM: reactants -> products
// with a bounds-aware arrow, the step's own metabolites emphasised, currency
// co-factors de-emphasised. The plain string (gemEquation) rides along as the
// accessible name.
function gemEquationHTML(g, rid, fromMid, toMid) {
  const r = g.reactions.find(x => x.id === rid);
  if (!r) return '';
  const term = (m, c) => {
    const coef = Math.abs(c) === 1 ? '' : `<span class="eq-coef">${Math.round(Math.abs(c) * 1e4) / 1e4}</span> `;
    const cur = !!(GRAPH && GRAPH.metabolites[m] && GRAPH.metabolites[m].cur);
    const main = m === fromMid || m === toMid;
    const nm = metName(m);
    return `<span class="eq-term${main ? ' eq-main' : cur ? ' eq-cur' : ''}" title="${esc(nm ? `${nm} (${m})` : m)}">${coef}${esc(m)}</span>`;
  };
  const lhs = [], rhs = [];
  for (const [m, c] of Object.entries(r.stoich)) (c < 0 ? lhs : rhs).push(term(m, c));
  const arrow = (r.lb < 0 && r.ub > 0) ? '⇌' : (r.ub <= 0 && r.lb < 0) ? '←' : '→';
  const plain = gemEquation(g, rid);
  return `<span class="eq-side">${lhs.join('<span class="eq-plus" aria-hidden="true">+</span>')}</span>` +
    `<span class="eq-arrow" aria-label="${esc(plain || '')}">${arrow}</span>` +
    `<span class="eq-side">${rhs.join('<span class="eq-plus" aria-hidden="true">+</span>')}</span>`;
}

// Subsystem span of a pathway in the context GEM: which subsystems its
// carried reactions belong to, as chips coloured by pathway group.
function subsystemSummaryHTML(g, pw) {
  const counts = new Map();   // subsystem -> {n, group}
  let carried = 0, totalAlts = 0;
  for (const st of pw.steps) for (const alt of st.rxns) {
    totalAlts++;
    const r = g.reactions.find(x => x.id === alt.id);
    if (!r) continue;
    carried++;
    const key = r.subsystem || '(no subsystem recorded)';
    if (!counts.has(key)) counts.set(key, { n: 0, group: r.group });
    counts.get(key).n++;
  }
  if (!carried) return '';
  const chips = [...counts.entries()].sort((a, b) => b[1].n - a[1].n).map(([s, v]) =>
    `<span class="subchip" style="--gc:${GROUP_COLORS[v.group] || '#98948C'}">${esc(subsystemLabel(s))} · ${v.n}</span>`).join('');
  return `<div class="pw-subsys">
    <span class="subsys-k">Subsystems spanned (${counts.size}, over ${carried} of ${totalAlts} candidate reactions in <span class="mono">${esc(getContext().gem)}</span>):</span>
    ${chips}
  </div>`;
}

// Default header facts for the Simulate columns of a card: honest absence
// until a run covers the pathway. mode2 overwrites them; a cleared run
// restores them (onSearchUpdate 'cleared').
const SIM_FACT_DEFAULTS = {
  feas: '<span class="fact-na">not tested; run Simulate (stage 3)</span>',
  flux: '<span class="fact-na">not computed</span>',
  yield: '<span class="fact-na">not computed</span>',
};

function resetSimFacts() {
  document.querySelectorAll('#results-body .pcard .pwh-facts [data-pwh]').forEach(el => {
    const dd = el.querySelector('dd');
    if (dd) dd.innerHTML = SIM_FACT_DEFAULTS[el.dataset.pwh] || '<span class="fact-na">not computed</span>';
  });
}

// One metabolite node on the flow rail.
function flowNodeHTML(mid, isEnd) {
  const nm = metName(mid);
  return `<div class="flow-node${isEnd ? ' flow-end' : ''}">
    <span class="fn-dot" aria-hidden="true"></span>
    <span class="fn-name">${esc(nm || mid)}</span>
    ${nm ? `<span class="fn-id mono">${esc(mid)}</span>` : ''}
  </div>`;
}

function pathwayCard(pw, idx) {
  const nGems = META.accs.length;
  const d = document.createElement('details');
  d.className = 'pcard';
  d.dataset.pwIdx = idx;
  const chainMids = pw.mets;
  const chainHtml = chainMids.map((m, i) => {
    const nm = metName(m) || m;
    const cls = i === 0 || i === chainMids.length - 1 ? 'chain-node chain-end' : 'chain-node chain-mid';
    return `<span class="${cls}">${esc(nm)}</span>`;
  }).join('<span class="chain-arrow" aria-hidden="true">→</span>');
  d.innerHTML = `
    <summary>
      <span class="rank">#${idx + 1}</span>
      <span class="plen">${pw.len} step${pw.len > 1 ? 's' : ''}</span>
      <span class="carriers"><span class="carrier-meter" aria-hidden="true"><span class="carrier-fill" style="width:${Math.round(pw.carriers / nGems * 100)}%"></span></span>${pw.carriers} of ${nGems} GEMs</span>
      <span class="feas-slot"></span>
      <span class="pcard-marker" aria-hidden="true"></span>
      <span class="chain" title="${esc(pw.mets.join(' → '))}">${chainHtml}</span>
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
  d.addEventListener('toggle', () => {
    if (d.open && map) showPathwayOnMap(pw, idx, { animate: true });
    if (d.open) enrichCard(d, pw);
  });

  // ---- header: route pills, facts with denominators, pathway-group key
  const header = document.createElement('header');
  header.className = 'pw-header';
  const routeHtml = chainMids.map((m, i) => {
    const nm = metName(m);
    const end = i === 0 || i === chainMids.length - 1;
    return `<span class="route-pill${end ? ' route-end' : ''}" title="${esc(m)}">${esc(nm || m)}</span>`;
  }).join('<span class="route-sep" aria-hidden="true">→</span>');
  const groupCounts = new Map();     // union pathway group -> candidate reaction count
  let nCand = 0;
  for (const st of pw.steps) for (const r of st.rxns) {
    nCand++;
    groupCounts.set(r.group, (groupCounts.get(r.group) || 0) + 1);
  }
  const groupKey = [...groupCounts.entries()].sort((a, b) => b[1] - a[1]).map(([g, n]) =>
    `<span class="gkey"><span class="gdot" style="background:${GROUP_COLORS[g] || '#98948C'}"></span>${esc(g)} · ${n}</span>`).join('');
  header.innerHTML = `
    <div class="pwh-route">${routeHtml}</div>
    <dl class="pwh-facts">
      <div class="pwh-fact"><dt>Length</dt><dd>${pw.len} step${pw.len > 1 ? 's' : ''}</dd></div>
      <div class="pwh-fact"><dt>Carried end-to-end</dt><dd>${pw.carriers} of ${nGems} GEMs</dd></div>
      <div class="pwh-fact" data-pwh="feas"><dt>Feasibility</dt><dd>${SIM_FACT_DEFAULTS.feas}</dd></div>
      <div class="pwh-fact" data-pwh="flux"><dt>Max product flux</dt><dd>${SIM_FACT_DEFAULTS.flux}</dd></div>
      <div class="pwh-fact" data-pwh="yield"><dt>Yield</dt><dd>${SIM_FACT_DEFAULTS.yield}</dd></div>
    </dl>
    <div class="pwh-keys">
      <span class="pwh-keyk">Pathway groups (${groupCounts.size}, over ${nCand} candidate reaction${nCand > 1 ? 's' : ''}):</span>
      ${groupKey}
    </div>
    <div class="pw-geminfo"><p class="termination">Choose a GEM (Model or Simulate stage) to see each reaction's genes, equation and subsystem in that strain.</p></div>`;
  bodyEl.appendChild(header);

  // ---- step flow: metabolite nodes on a rail, one designed card per step
  const flow = document.createElement('div');
  flow.className = 'pw-flow';
  const altMax = 6;
  let flowHtml = '';
  pw.steps.forEach((st, i) => {
    flowHtml += flowNodeHTML(st.from, i === 0);
    const shown = st.rxns.slice(0, altMax);
    // step accent: the group colour when every shown candidate agrees, else neutral
    const stepGroups = new Set(shown.map(r => r.group));
    const sg = stepGroups.size === 1 ? (GROUP_COLORS[shown[0].group] || 'var(--line-strong)') : 'var(--line-strong)';
    const rows = shown.map(r => `
      <div class="fs-rxn rxn-alt" data-rxn="${esc(r.id)}" data-from="${esc(st.from)}" data-to="${esc(st.to)}">
        <div class="fs-info">
          <div class="fs-enzyme"><span class="gs-genes"></span></div>
          <div class="fs-idline">
            <button class="rid-btn" type="button" data-rid="${esc(r.id)}" title="Open ${esc(r.id)} in the Model stage">${esc(r.id)}</button>
            ${r.name && r.name !== r.id ? `<span class="fs-rname">${esc(r.name)}</span>` : ''}
          </div>
          <div class="fs-tags">
            <span class="subchip" style="--gc:${GROUP_COLORS[r.group] || '#98948C'}">${esc(r.group)}</span>
            <span class="gs-sub"></span>
            ${r.ec.length ? `<span class="fs-ec">EC ${esc(r.ec.join(', '))}</span>` : ''}
            ${r.dir === 'rev' ? '<span class="fs-dir">runs reverse of written direction</span>' : ''}
          </div>
          <div class="fs-eq gs-eq"></div>
        </div>
        <div class="fs-flux"><div class="flux-slot"></div></div>
      </div>`).join('');
    flowHtml += `
      <article class="flow-step" data-step="${i}" style="--sg:${sg}" aria-label="Step ${i + 1} of ${pw.len}: ${esc(metName(st.from) || st.from)} to ${esc(metName(st.to) || st.to)}">
        <div class="fs-head">
          <span class="fs-no">Step ${i + 1} of ${pw.len}</span>
          ${st.rxns.length > 1 ? `<span class="fs-altnote">any 1 of ${st.rxns.length} candidate reactions carries this step</span>` : ''}
        </div>
        <div class="fs-rxns">${rows}</div>
        ${st.rxns.length > altMax ? `<p class="termination">Showing ${altMax} of ${st.rxns.length} candidate reactions for this step.</p>` : ''}
      </article>`;
  });
  flowHtml += flowNodeHTML(pw.mets[pw.mets.length - 1], true);
  flow.innerHTML = flowHtml;
  bodyEl.appendChild(flow);

  // ---- side column: strain presence per species
  const side = document.createElement('aside');
  side.className = 'pw-side';
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
  const sideH = document.createElement('h3');
  sideH.className = 'side-h';
  sideH.textContent = `Strain presence (${pw.carriers} of ${nGems} GEMs carry every step)`;
  const note = document.createElement('p');
  note.className = 'termination';
  note.textContent = 'A filled cell means the strain’s GEM contains at least one reaction for every step. Whether the strain’s own bounds and a medium permit flux through every step is a Simulate check.';
  side.append(sideH, matrix, note);
  bodyEl.appendChild(side);

  // Feasibility detail slot (filled by the Simulate stage when a run covers this card)
  const m2slot = document.createElement('div');
  m2slot.className = 'mode2-slot';
  bodyEl.appendChild(m2slot);

  // actions
  const actions = document.createElement('div');
  actions.className = 'cardactions';
  const actK = document.createElement('span');
  actK.className = 'actions-k';
  actK.textContent = 'This pathway:';
  actions.appendChild(actK);
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

  bodyEl.querySelectorAll('button.rid-btn').forEach(b => b.addEventListener('click', () => {
    openReactionInBrowser(b.dataset.rid);
    location.hash = '#/model';
  }));
  return d;
}

// Cards already open when a context GEM lands (chosen in Model or Simulate
// after expanding) get their gene / equation / subsystem slots filled without
// needing a close-and-reopen.
function enrichOpenCards() {
  if (!lastResults) return;
  document.querySelectorAll('#results-body .pcard[open]').forEach(d => {
    const pw = lastResults.res.pathways[+d.dataset.pwIdx];
    if (pw) enrichCard(d, pw);
  });
}

// Fill a card's per-reaction gene / equation / subsystem slots from the
// context GEM. Runs on open and re-runs when the context GEM changes; without
// a chosen GEM the prompt stays and nothing is invented.
async function enrichCard(d, pw) {
  const acc = getContext().gem;
  if (!acc) return;
  if (d.dataset.enrichedAcc === acc) return;
  if (!ctxGem || ctxGem.acc !== acc) {
    try { ctxGem = await loadGem(acc); } catch { return; }
    if (getContext().gem !== acc) return;
  }
  d.dataset.enrichedAcc = acc;
  d.querySelectorAll('.rxn-alt').forEach(row => {
    const rid = row.dataset.rxn;
    const genesEl = row.querySelector('.gs-genes');
    const subEl = row.querySelector('.gs-sub');
    const eqEl = row.querySelector('.gs-eq');
    if (!genesEl || !subEl || !eqEl) return;
    const r = ctxGem.reactions.find(x => x.id === rid);
    if (!r) {
      row.classList.add('rxn-absent');
      genesEl.innerHTML = '';
      subEl.innerHTML = `<span class="absent-tag">not in ${esc(acc)}</span>`;
      eqEl.innerHTML = '';
      return;
    }
    row.classList.remove('rxn-absent');
    const genes = (r.genes || []);
    genesEl.innerHTML = genes.length
      ? genes.map(l => geneLabelHTML(ctxGem, l)).join('<span class="gsep" aria-hidden="true">·</span>')
      : `<span class="no-gene">no gene rule in ${esc(acc)}</span>`;
    // the GEM's raw subsystem, unless it only repeats the union group chip
    const subLabel = r.subsystem ? subsystemLabel(r.subsystem) : '';
    const groupChipText = row.querySelector('.fs-tags .subchip');
    subEl.innerHTML = (subLabel && !(groupChipText && groupChipText.textContent.trim() === subLabel))
      ? `<span class="subchip" style="--gc:${GROUP_COLORS[r.group] || '#98948C'}">${esc(subLabel)}</span>` : '';
    eqEl.innerHTML = gemEquationHTML(ctxGem, rid, row.dataset.from, row.dataset.to);
  });
  const gemInfo = d.querySelector('.pw-geminfo');
  if (gemInfo) gemInfo.innerHTML = subsystemSummaryHTML(ctxGem, pw) ||
    `<p class="termination"><span class="mono">${esc(acc)}</span> carries none of this pathway's candidate reactions; genes and subsystems are shown per carried reaction only.</p>`;
}

function exportPathway(pw, idx, kind) {
  const presence = META.accs.map((a, i) => ({
    acc: a.acc, species: a.sp, carries_every_step: ((pw.mask >> BigInt(i)) & 1n) === 1n,
  }));
  const run = runView();
  const r = run ? run.results.get(idx) : null;
  const deep = run ? run.deep.get(idx) : null;
  const acc = getContext().gem;
  const g = (ctxGem && ctxGem.acc === acc) ? ctxGem : null;
  // genes as "symbol (locus)" for the context GEM; absent stays absent
  const genesOf = (rid) => {
    if (!g) return null;
    const gr = g.reactions.find(x => x.id === rid);
    return gr ? (gr.genes || []).map(l => geneLabel(g, l)) : null;
  };
  const meta = {
    substrate: lastResults.sub, product: lastResults.prod,
    species_filter: lastResults.species, rank: idx + 1, steps: pw.len,
    carriers: `${pw.carriers} of ${META.accs.length} GEMs`,
    gene_labels: g ? `symbol (locus) for GEM ${acc}` : 'not included (no GEM chosen)',
    flux_results: r
      ? `feasibility${deep ? ', pFBA, FVA and sampling' : ''} for ${run.acc} on ${run.medium}`
      : 'not included (no Simulate run covers this pathway)',
  };
  if (kind === 'json') {
    downloadBlob(JSON.stringify({
      ...meta,
      metabolites: pw.mets,
      reactions_per_step: pw.steps.map(s => ({
        from: s.from, to: s.to,
        reactions: s.rxns.map(x => ({
          id: x.id, name: x.name, direction: x.dir, ec: x.ec,
          genes: genesOf(x.id) ?? undefined,
        })),
      })),
      presence,
      flux: r ? {
        gem: run.acc, medium: run.medium,
        testable: r.testable, feasible: r.testable ? r.feasible : null,
        max_product_flux_mmol_gDW_h: r.testable && r.feasible ? r.productFlux : null,
        fva_ranges: deep && deep.fva ? deep.fva.ranges : null,
        sampling: deep && deep.sample && deep.sample.ok ? {
          samples: deep.sample.samples, requested: deep.stamp.requested, failed: deep.sample.failed,
          seed: deep.sample.seed, sampler: deep.sample.sampler,
          stats: Object.fromEntries(deep.rids.map(rid => [rid, deep.sample.stats.get(rid) || null])),
        } : null,
      } : null,
    }, null, 1), `pathway_${idx + 1}_${lastResults.sub}_to_${lastResults.prod}.json`, 'application/json');
  } else {
    const lines = [
      `# pathway ${idx + 1}: ${lastResults.sub} to ${lastResults.prod}; ${pw.len} steps; carried by ${pw.carriers} of ${META.accs.length} GEMs`,
      `# genes: ${meta.gene_labels}; flux: ${meta.flux_results}`,
      'step,from,to,reaction_ids,directions,ec_numbers,genes',
      ...pw.steps.map((s, i) => [i + 1, s.from, s.to,
        s.rxns.map(x => x.id).join(';'), s.rxns.map(x => x.dir).join(';'),
        s.rxns.flatMap(x => x.ec).join(';'),
        s.rxns.map(x => (genesOf(x.id) || []).join('|')).join(';'),
      ].map(csvEscape).join(',')),
      '',
      'acc,species,carries_every_step',
      ...presence.map(p => [p.acc, p.species, p.carries_every_step].map(csvEscape).join(',')),
    ];
    if (deep && deep.fva) {
      lines.push('', 'reaction,fva_min,fva_max,sample_median,sample_p5,sample_p95,sample_min,sample_max');
      for (const rid of deep.rids) {
        const rr = deep.fva.ranges[rid] || {};
        const ss = (deep.sample && deep.sample.ok && deep.sample.stats.get(rid)) || {};
        lines.push([rid, rr.min, rr.max, ss.median, ss.p5, ss.p95, ss.min, ss.max].map(csvEscape).join(','));
      }
    }
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
  // (the representative biomass assembly is a map visual and is not counted)
  const groups = GRAPH.groups.filter(g => g !== 'Biomass');
  const spNames = INDEX.species.map(s => s.name);
  const counts = spNames.map(() => groups.map(() => 0));
  const groupTotals = groups.map(() => 0);
  for (const r of GRAPH.reactions) {
    if (r.core) continue;
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
// With a feasibility run active the operator's colour tiers apply: the
// selected pathway is vivid cyan when feasible (species accent when not),
// every other feasible pathway is pale cyan, infeasible candidates stay dim.
function showPathwayOnMap(pw, idx, { animate = false, weights = null } = {}) {
  if (!map) return;
  highlightedIdx = idx;
  const run = runView();
  if (run && run.feasibleIdx.length) {
    const isFeasible = run.feasibleIdx.includes(idx);
    const others = run.feasibleIdx.filter(i => i !== idx)
      .map(i => lastResults.res.pathways[i].mets);
    let w = weights;
    if (!w && idx === run.bestIdx && run.bestWeights) w = run.bestWeights;
    map.highlightPathway(pw.mets, isFeasible ? SEARCH_COLORS.best : accentInk(), {
      animate, weights: w, others, altColor: SEARCH_COLORS.alt,
    });
  } else {
    map.highlightPathway(pw.mets, accentInk(), { animate, weights });
  }
  $('#clear-highlight').hidden = false;
}

// A feasibility / FVA / sampling run landed (or was cleared): refresh the
// length histogram, redraw the tiers, and draw the flux-carrying layer from
// the best pathway's pFBA solution.
function onSearchUpdate(phase) {
  if (phase === 'cleared') {
    if (map) { map.clearFluxEdges(); }
    fluxLayerOn = false;
    if (lastResults) { renderListControls(); }
    resetSimFacts();
    updateSearchLegend(null);
    return;
  }
  if (!lastResults) return;
  const run = runView();
  if (phase === 'done') { renderListControls(); }
  if (!map || !run) return;
  let fluxDrawn = null;
  if (run.fluxRids && !fluxLayerOn) {
    fluxDrawn = map.setFluxEdges(run.fluxRids, SEARCH_COLORS.fluxEdge);
    fluxLayerOn = true;
  }
  if (run.feasibleIdx.length) {
    const cur = (highlightedIdx != null && run.feasibleIdx.includes(highlightedIdx))
      ? highlightedIdx : run.bestIdx;
    showPathwayOnMap(lastResults.res.pathways[cur], cur, { animate: false });
  }
  updateSearchLegend(run, fluxDrawn);
}

// Search-tier legend block inside the map legend, present only while the
// tiers themselves are drawn (a run with zero feasible pathways draws none).
let searchLegendData = null;
function updateSearchLegend(run, fluxDrawn) {
  searchLegendData = (run && (run.feasibleIdx.length || run.fluxRids))
    ? { nFeas: run.feasibleIdx.length, flux: !!run.fluxRids, fluxDrawn } : null;
  buildLegend();
}

function clearMapHighlight() {
  if (map) { map.clearHighlight(); map.clearReactionHighlight(); map.clearFluxEdges(); }
  fluxLayerOn = false;
  highlightedIdx = null;
  updateSearchLegend(null);
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
      resultsSummary: () => $('#results-summary'),
      speciesTokens: SPECIES_TOKENS,
      onStateChange: () => refreshSimUI(),
      onSearchUpdate,
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
    to ${metLabelHTML(lastResults.prod)} ready: feasibility (LP) on the shortest ${Math.min(10, n)} of ${fmt.format(n)},
    then FVA + flux sampling on the ${Math.min(5, n)} shortest feasible; the rest per card.`;
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
  if (!el || !GRAPH) return;
  el.hidden = false;
  const searchPart = searchLegendData ? (
    '<span class="lg" style="font-weight:600">Flux search:</span>' +
    (searchLegendData.nFeas > 0 ? `<span class="lg"><span class="swatch" style="background:${SEARCH_COLORS.best}"></span>selected / shortest feasible pathway</span>` : '') +
    (searchLegendData.nFeas > 1 ? `<span class="lg"><span class="swatch" style="background:${SEARCH_COLORS.alt}"></span>other feasible pathways (${searchLegendData.nFeas - 1})</span>` : '') +
    (searchLegendData.flux ? `<span class="lg"><span class="swatch" style="background:${SEARCH_COLORS.fluxEdge}"></span>carries flux in the best pathway's pFBA optimum (currency edges not drawn)</span>` : '') +
    '<span class="lg">dim = infeasible or untested</span>'
  ) : '';
  const nEsch = Object.values(GRAPH.metabolites).filter(m => m.esch).length;
  el.innerHTML =
    searchPart +
    '<span class="lg" style="font-weight:600">Shells:</span>' +
    '<span class="lg">outer = extracellular + exchange</span>' +
    '<span class="lg">middle = periplasm</span>' +
    '<span class="lg">inner = cytosol</span>' +
    (nEsch ? `<span class="lg">flat central sheet = ${nEsch} of ${fmt.format(GRAPH.n_metabolites)} metabolites (central carbon) at the Escher e_coli_core layout</span>` : '') +
    '<span class="lg" style="font-weight:600;margin-left:6px">Groups:</span>' + GRAPH.groups.map(g =>
      `<span class="lg"><span class="swatch" style="background:${GROUP_COLORS[g] || '#98948C'}"></span>${esc(g)}</span>`).join('') +
    '<span class="lg"><span class="swatch" style="background:' + (GROUP_COLORS.Biomass || '#3A3E45') + ';border-radius:50%"></span>central node = representative biomass (93 common precursors)</span>';
}
