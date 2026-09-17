// Metabolic Atlas: main controller. Loads the registry + union graph, wires the
// substrate/product pickers, species filter, Mode-1 search, the 3D map and the
// GEM browser. Absent values render as absent; every count carries its denominator.

import { loadIndex, loadGraph, loadGraphMeta, fmt, downloadBlob, csvEscape } from './data.js';
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
let metList = [];                        // [{mid, name, comp, group, cur}]
let map = null;                          // map3d api or null
let selectedSpecies = new Set();
let lastResults = null;

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

  GROUP_COLORS = groupColorMap(GRAPH.groups);
  metList = Object.entries(GRAPH.metabolites).map(([mid, m]) => ({
    mid, name: m.n !== mid ? m.n : '', comp: m.c, group: m.g, cur: !!m.cur,
  }));
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

  function renderOptions(q) {
    const query = q.trim().toLowerCase();
    if (!query) { close(); return; }
    const starts = [], contains = [];
    for (const m of metList) {
      const idL = m.mid.toLowerCase(), nmL = m.name.toLowerCase();
      if (idL.startsWith(query) || nmL.startsWith(query)) starts.push(m);
      else if (idL.includes(query) || nmL.includes(query)) contains.push(m);
    }
    const all = starts.concat(contains);
    options = all.slice(0, 50);
    if (!all.length) {
      listbox.innerHTML = `<li class="mcap" role="presentation">No metabolite matches "${esc(q)}" among ${fmt.format(metList.length)} ids in the union map.</li>`;
      open(); return;
    }
    const cap = all.length > options.length
      ? `<li class="mcap" role="presentation">Showing ${options.length} of ${fmt.format(all.length)} matches; keep typing to narrow.</li>` : '';
    listbox.innerHTML = options.map((m, i) => `
      <li id="${kind}-opt-${i}" role="option" aria-selected="false" data-mid="${esc(m.mid)}">
        <span class="mid">${esc(m.mid)}</span>
        ${m.name ? `<span class="mname">${esc(m.name)}</span>` : ''}
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
    input.value = mid;
    close();
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
      pickerState.sub.mid = a; $('#sub-input').value = a;
      pickerState.prod.mid = b; $('#prod-input').value = b;
      runSearch();
    });
    row.appendChild(btn);
  }
}

// ================= Mode-1 search =================
function runSearch() {
  const sub = pickerState.sub.mid, prod = pickerState.prod.mid;
  const summary = $('#results-summary');
  if (!GRAPH || !META) return;
  if (!sub || !prod) {
    summary.innerHTML = `<span class="status">Pick both a substrate and a product from the list (type to search ${fmt.format(metList.length)} metabolite ids).</span>`;
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
      <span class="mono">${esc(sub)}</span> to <span class="mono">${esc(prod)}</span>
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
    <span class="mono">${esc(sub)}</span> to <span class="mono">${esc(prod)}</span> ·
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
  const chain = pw.mets.join(' → ');
  d.innerHTML = `
    <summary>
      <span class="rank">#${idx + 1}</span>
      <span class="plen">${pw.len} step${pw.len > 1 ? 's' : ''}</span>
      <span class="carriers">${pw.carriers} of ${nGems} GEMs carry every step</span>
      <span class="chain">${esc(chain)}</span>
    </summary>
    <div class="pcard-body"></div>`;
  const bodyEl = d.querySelector('.pcard-body');

  const stepsEl = document.createElement('div');
  pw.steps.forEach((st, i) => {
    const CAP = 6;
    const shown = st.rxns.slice(0, CAP);
    const stepDiv = document.createElement('div');
    stepDiv.className = 'step';
    stepDiv.innerHTML = `<div class="step-mets">Step ${i + 1}: ${esc(st.from)} → ${esc(st.to)}</div>` +
      shown.map(r => `
        <div class="rxn-alt">
          <button class="rid btn small" type="button" data-rid="${esc(r.id)}" title="Open in GEM browser">${esc(r.id)}</button>
          ${r.dir === 'rev' ? '<span class="ec">(reverse of written direction)</span>' : ''}
          ${r.ec.length ? `<span class="ec">EC ${esc(r.ec.join(', '))}</span>` : ''}
          <span class="rn">${esc(r.name || '')}</span>
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
  note.textContent = 'A filled cell means the strain’s GEM contains at least one reaction for every step. Per-strain direction feasibility is a Mode 2 (phase 2) check.';
  matrix.appendChild(note);
  bodyEl.appendChild(matrix);

  // actions
  const actions = document.createElement('div');
  actions.className = 'cardactions';
  const btnMap = document.createElement('button');
  btnMap.className = 'btn small'; btnMap.type = 'button';
  btnMap.textContent = map ? 'Show on 3D map' : '3D map unavailable';
  btnMap.disabled = !map;
  btnMap.addEventListener('click', () => {
    if (!map) return;
    const ink = getComputedStyle(document.documentElement).getPropertyValue('--accent-ink').trim() || '#245C77';
    map.highlightPathway(pw.mets, ink);
    $('#clear-highlight').hidden = false;
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
    note: 'Mode 1 enzyme presence; per-strain direction feasibility not checked (phase 2).',
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
        tip.innerHTML = `<span class="tid">${esc(hit.mid)}</span><br>
          ${hit.met.n !== hit.mid ? esc(hit.met.n) + '<br>' : ''}
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
    $('#clear-highlight').addEventListener('click', () => { map.clearHighlight(); $('#clear-highlight').hidden = true; });
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
