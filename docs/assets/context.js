// Session context: the four choices that carry across every stage.
//   sub / prod  - substrate and product metabolite ids (set in Discover)
//   gem         - GEM accession (set in Model, Simulate or Engineer)
//   medium      - medium label (set in Simulate or Engineer)
// One store, one writer function, subscribers per stage. An unset value is
// null and renders as "not set", never as a default the user did not choose.

const state = { sub: null, prod: null, gem: null, medium: null };
const subs = new Set();

export function getContext() { return { ...state }; }

export function setContext(patch) {
  const changed = [];
  for (const k of ['sub', 'prod', 'gem', 'medium']) {
    if (k in patch && patch[k] !== state[k]) { state[k] = patch[k]; changed.push(k); }
  }
  if (changed.length) for (const cb of subs) cb(getContext(), changed);
}

export function onContext(cb) { subs.add(cb); return () => subs.delete(cb); }

// ---- context bar: always visible under the header; each item links to the
// stage that sets it.
let barEl = null, labelers = null;

export function initContextBar(el, opts) {
  barEl = el;
  labelers = opts;   // { metLabel(mid), gemLabel(acc) }
  onContext(renderBar);
  renderBar();
}

const esc = (s) => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

function item(href, key, value) {
  return `<a class="ctx-item ${value ? 'set' : 'unset'}" href="${href}">
    <span class="ctx-k">${key}</span>
    <span class="ctx-v">${value ? esc(value) : 'not set'}</span></a>`;
}

function renderBar() {
  if (!barEl) return;
  const met = (mid) => (labelers && labelers.metLabel(mid)) || mid;
  const gem = (acc) => (labelers && labelers.gemLabel(acc)) || acc;
  barEl.innerHTML =
    item('#/discover', 'Substrate', state.sub ? met(state.sub) : null) +
    '<span class="ctx-arrow" aria-hidden="true">&#8594;</span>' +
    item('#/discover', 'Product', state.prod ? met(state.prod) : null) +
    item('#/model', 'GEM', state.gem ? gem(state.gem) : null) +
    item('#/simulate', 'Medium', state.medium);
}
