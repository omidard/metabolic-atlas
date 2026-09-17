// Metabolic Atlas: main controller. Loads the registry + union graph, wires the
// substrate/product pickers, species filter, Mode-1 search, the 3D map and the
// GEM browser. Absent values render as absent; every count carries its denominator.

import { loadIndex, loadGraph, loadGraphMeta, loadMetIndex, fmt, downloadBlob, csvEscape } from './data.js';
import { searchPathways, carriersBySpecies, DEFAULTS } from './search.js';
import { initGemView, openReactionInBrowser, setIndexData } from './gem.js';

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

// ---- pathway-group colours: 13 distinct hues + neutral for Other
function groupColorMap(groups) {
  const named = groups.filter(g => g !== 'Other');
  const out = {};
  named.forEach((g, i) => {
    const h = Math.round(i * 360 / named.length);
    out[g] = hslToHex(h, 52, 40);
  });
  out['Other'] = '#98948C';
  return out;
}
function hslToHex(h, s, l) {
  s /= 100; l /= 100;
  const k = n => (n + h / 30) % 12;
  const a = s * Math.min(l, 1 - l);
  const f = n => l - a * Math.max(-1, Math.min(k(n) - 3, Math.min(9 - k(n), 1)));
  return '#' + [f(0), f(8), f(4)].map(x => Math.round(255 * x).toString(16).padStart(2, '0')).join('');
}

// ---- app state
let INDEX = null, GRAPH = null, META = null, GROUP_COLORS = null;
let MET_INDEX = null;                    // {mid: {name, kegg, chebi}} or null when unavailable
let metList = [];                        // [{mid, name, kegg, comp, group, cur}]
let map = null;                          // map3d api or null
let selectedSpecies = new Set();
let lastResults = null;
let MODE2 = null;                        // lazily imported mode2 module api or null
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

  try {
    INDEX = await loadIndex();
  } catch (e) {
    $('#dataset-line').textContent = `Could not load the GEM index (${e.message}). Reload the page to retry.`;
    return;
  }
  setIndexData(INDEX);
  initGemView($('#view-gems'), INDEX, { onAccentChange: setAccent });
  buildSpeciesChips();

  try {
    [GRAPH, META] = await Promise.all([loadGraph(), loadGraphMeta()]);
  } catch (e) {
    $('#dataset-line').textContent = `Could not load the union map (${e.message}). Reload the page to retry; the GEM browser works without it.`;
    $('#map-status').textContent = 'Union map unavailable. Reload the page to retry.';
    return;
  }

  try {
    MET_INDEX = await loadMetIndex();
  } catch (e) {
    MET_INDEX = null;   // pickers fall back to id-only matching and say so
  }

  GROUP_COLORS = groupColorMap(GRAPH.groups);
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
  wireModeToggle();

  initMap();
})();

// ================= routing =================
function route() {
  const h = location.hash || '#/search';
  const onGems = h.startsWith('#/gems');
  $('#view-search').hidden = onGems;
  $('#view-gems').hidden = !onGems;
  $('#nav-search').setAttribute('aria-current', onGems ? 'false' : 'page');
  $('#nav-gems').setAttribute('aria-current', onGems ? 'page' : 'false');
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

// ================= search (Mode 1 enumeration; Mode 2 adds feasibility) =================
function currentMode() {
  const r = document.querySelector('input[name="mode"]:checked');
  return r ? r.value : '1';
}

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
  const mode = currentMode();
  if (mode === '2' && MODE2) {
    const why = MODE2.notReadyReason();
    if (why) {
      summary.innerHTML = `<span class="status">${esc(why)}</span>`;
      return;
    }
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
    if (mode === '2' && MODE2) MODE2.runFeasibility(lastResults);
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
    const CAP = 6;
    const shown = st.rxns.slice(0, CAP);
    const stepDiv = document.createElement('div');
    stepDiv.className = 'step';
    stepDiv.dataset.step = i;
    const fromNm = metName(st.from), toNm = metName(st.to);
    stepDiv.innerHTML = `<div class="step-mets">Step ${i + 1}: ${fromNm ? esc(fromNm) + ' ' : ''}<span class="mono">${esc(st.from)}</span> → ${toNm ? esc(toNm) + ' ' : ''}<span class="mono">${esc(st.to)}</span></div>` +
      shown.map(r => `
        <div class="rxn-alt" data-rxn="${esc(r.id)}">
          <button class="rid btn small" type="button" data-rid="${esc(r.id)}" title="Open in GEM browser">${esc(r.id)}</button>
          ${r.dir === 'rev' ? '<span class="ec">(reverse of written direction)</span>' : ''}
          ${r.ec.length ? `<span class="ec">EC ${esc(r.ec.join(', '))}</span>` : ''}
          <span class="rn">${esc(r.name || '')}</span>
          <span class="flux-slot"></span>
        </div>`).join('') +
      (st.rxns.length > CAP ? `<p class="termination">Showing ${CAP} of ${st.rxns.length} alternative reactions for this step.</p>` : '');
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
  note.textContent = 'A filled cell means the strain’s GEM contains at least one reaction for every step. Whether the strain’s own bounds and a medium permit flux through every step is a Mode 2 check.';
  matrix.appendChild(note);
  bodyEl.appendChild(matrix);

  // Mode-2 feasibility detail slot (filled by mode2.js when a run covers this card)
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
    location.hash = '#/gems';
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
    note: 'Mode 1 enzyme presence; flux feasibility under a medium is a Mode 2 check and is not included in this export.',
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
  if (map) map.clearHighlight();
  highlightedIdx = null;
  $('#clear-highlight').hidden = true;
}

// ================= Mode 2 module (lazy) =================
async function activateMode2() {
  if (MODE2) { MODE2.setActive(true); return; }
  const panel = $('#mode2-panel');
  panel.hidden = false;
  panel.innerHTML = '<div class="card"><span class="status">Loading the flux engine (GLPK-WASM, about 250 kB)…</span></div>';
  try {
    const mod = await import('./mode2.js');
    MODE2 = await mod.initMode2(panel, {
      index: INDEX, graphMeta: META,
      metName, metLabelHTML, setAccent,
      getEndpoints: () => ({ sub: pickerState.sub.mid, prod: pickerState.prod.mid }),
      findCard: (i) => document.querySelector(`.pcard[data-pw-idx="${i}"]`),
      showFluxOnMap: (pw, i, weights) => showPathwayOnMap(pw, i, { animate: false, weights }),
      resultsSummary: () => $('#results-summary'),
      speciesTokens: SPECIES_TOKENS,
    });
    MODE2.onEndpointsChange(pickerState.sub.mid, pickerState.prod.mid);
    MODE2.setActive(true);
  } catch (e) {
    panel.innerHTML = `<div class="card"><p class="status error">Mode 2 is unavailable: the flux engine failed to load (${esc(e.message)}).
      Reload the page to retry. Mode 1 search keeps working.</p></div>`;
    const m1 = document.querySelector('input[name="mode"][value="1"]');
    if (m1) m1.checked = true;
  }
}

function wireModeToggle() {
  document.querySelectorAll('input[name="mode"]').forEach(r => {
    r.addEventListener('change', () => {
      const m2 = currentMode() === '2';
      $('#mode2-panel').hidden = !m2;
      if (m2) activateMode2();
      else if (MODE2) MODE2.setActive(false);
    });
  });
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
      Pathway search and the GEM browser work without the 3D view.`;
  }
}

function buildLegend() {
  const el = $('#map-legend');
  el.hidden = false;
  el.innerHTML = '<span class="lg" style="font-weight:600">Pathway groups:</span>' + GRAPH.groups.map(g =>
    `<span class="lg"><span class="swatch" style="background:${GROUP_COLORS[g]}"></span>${esc(g)}</span>`).join('');
}
