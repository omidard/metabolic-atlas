// Data loading and caching. All fetches are relative to docs/ root.
// Absent values stay null; loaders never invent a 0.

const cache = new Map();

async function fetchJSON(url, label, onStatus) {
  if (cache.has(url)) return cache.get(url);
  if (onStatus) onStatus(`Loading ${label}…`);
  const res = await fetch(url);
  if (!res.ok) {
    throw new Error(`${label}: HTTP ${res.status} loading ${url}`);
  }
  const data = await res.json();
  cache.set(url, data);
  return data;
}

export function loadIndex(onStatus) {
  return fetchJSON('data/index.json', 'GEM index', onStatus);
}

export function loadGraph(onStatus) {
  return fetchJSON('data/global_graph.json', 'union map (1.3 MB)', onStatus);
}

export function loadGraphMeta(onStatus) {
  return fetchJSON('assets/graph_meta.json', 'reaction presence index', onStatus);
}

export function loadGem(acc, onStatus) {
  return fetchJSON(`data/gems/${acc}.json`, `GEM ${acc}`, onStatus);
}

export function loadSeqs(acc, onStatus) {
  return fetchJSON(`data/seqs/${acc}.json`, `sequences for ${acc}`, onStatus);
}

export function loadMedia(onStatus) {
  return fetchJSON('data/media.json', 'media definitions', onStatus);
}

export function loadMetIndex(onStatus) {
  return fetchJSON('data/met_index.json', 'metabolite name index', onStatus);
}

// ---- in-memory bound edits: acc -> rxnId -> {lb, ub} (Phase 2 analyses read these)
const boundEdits = new Map();

export function getEdit(acc, rxnId) {
  const m = boundEdits.get(acc);
  return m ? (m.get(rxnId) || null) : null;
}

export function setEdit(acc, rxnId, lb, ub) {
  if (!boundEdits.has(acc)) boundEdits.set(acc, new Map());
  boundEdits.get(acc).set(rxnId, { lb, ub });
}

export function clearEdit(acc, rxnId) {
  const m = boundEdits.get(acc);
  if (m) m.delete(rxnId);
}

export function clearAllEdits(acc) {
  boundEdits.delete(acc);
}

export function editCount(acc) {
  const m = boundEdits.get(acc);
  return m ? m.size : 0;
}

// ---- shared formatting helpers
export const fmt = new Intl.NumberFormat('en-US');

export function n(x) {
  // absent renders as absent, never as a confident 0
  return (x === null || x === undefined) ? 'not computed' : fmt.format(x);
}

export function downloadBlob(text, filename, type) {
  const blob = new Blob([text], { type });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 5000);
}

export function csvEscape(v) {
  if (v === null || v === undefined) return '';
  const s = String(v);
  return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
}
