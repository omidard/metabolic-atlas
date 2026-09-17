// Mode 1: bioconversion pathway enumeration over the union reaction graph.
// Metabolite-level directed graph; each edge carries the set of reactions that
// realise it. Currency metabolites (cur:1) are never used as intermediate hubs;
// they may be a chosen source or target. Reaction direction comes from the union
// of per-GEM bounds (graph_meta: d=1 forward as written, 2 reverse only, 3 both).

export const DEFAULTS = { maxDepth: 8, maxPaths: 60, maxExpansions: 500000 };

// Build adjacency restricted to reactions carried by at least one selected species.
export function buildSearchGraph(graph, meta, selectedSpecies) {
  const spSet = new Set(selectedSpecies);
  const adj = new Map();        // met -> Map(nextMet -> rxn list)
  const addEdge = (u, v, rxn) => {
    if (u === v) return;
    let m = adj.get(u);
    if (!m) { m = new Map(); adj.set(u, m); }
    let list = m.get(v);
    if (!list) { list = []; m.set(v, list); }
    list.push(rxn);
  };
  for (const r of graph.reactions) {
    if (!r.sp.some(s => spSet.has(s))) continue;
    const info = meta.rxns[r.id] || null;
    const d = info ? info.d : 1;
    const rec = { id: r.id, name: r.n, group: r.g, ec: (info && info.e) || [], mask: info ? BigInt('0x' + info.m) : 0n };
    if (d & 1) for (const s of r.s) for (const p of r.p) addEdge(s, p, { ...rec, dir: 'fwd' });
    if (d & 2) for (const s of r.s) for (const p of r.p) addEdge(p, s, { ...rec, dir: 'rev' });
  }
  return adj;
}

export function searchPathways(graph, meta, source, target, selectedSpecies, opts = {}) {
  const { maxDepth, maxPaths, maxExpansions } = { ...DEFAULTS, ...opts };
  const mets = graph.metabolites;
  // Compartment forms of the chosen endpoints stay traversable even when they
  // carry the currency flag (a CO2 search must be able to cross the membrane).
  const baseOf = (mid) => {
    const c = mets[mid] && mets[mid].c;
    return c && mid.endsWith('_' + c) ? mid.slice(0, -(c.length + 1)) : mid;
  };
  const endpointBases = new Set([baseOf(source), baseOf(target)]);
  const isCurrency = (m) => !!(mets[m] && mets[m].cur) && !endpointBases.has(baseOf(m));
  const adj = buildSearchGraph(graph, meta, selectedSpecies);

  // Reverse BFS from the target: dist[v] = fewest edges from v to target,
  // never relaying through a currency metabolite. Used to prune the forward
  // enumeration to nodes that can still reach the target within maxDepth.
  const dist = new Map([[target, 0]]);
  {
    const radj = new Map();      // v -> Set(u) for every edge u -> v
    for (const [u, m] of adj) for (const v of m.keys()) {
      let s = radj.get(v);
      if (!s) { s = new Set(); radj.set(v, s); }
      s.add(u);
    }
    const q = [target];
    let qi = 0;
    while (qi < q.length) {
      const v = q[qi++];
      const dv = dist.get(v);
      if (dv >= maxDepth) continue;
      if (v !== target && isCurrency(v)) continue;   // currency never relays
      for (const u of radj.get(v) || []) {
        if (!dist.has(u)) { dist.set(u, dv + 1); q.push(u); }
      }
    }
  }

  const found = [];              // arrays of met ids
  let expansions = 0;
  let hitPaths = false, hitBudget = false;

  // BFS over simple paths: paths surface in non-decreasing length, so a hit
  // path cap keeps the shortest pathways, never a depth-first arbitrary subset.
  const queue = [[source]];
  let head = 0;
  outer: while (head < queue.length) {
    const path = queue[head++];
    const u = path[path.length - 1];
    const depth = path.length - 1;               // edges so far
    const next = adj.get(u);
    if (!next) continue;
    for (const v of next.keys()) {
      expansions++;
      if (expansions >= maxExpansions) { hitBudget = true; break outer; }
      if (path.includes(v)) continue;            // simple paths only
      if (v === target) {
        if (depth + 1 <= maxDepth) {
          found.push([...path, v]);
          if (found.length >= maxPaths) { hitPaths = true; break outer; }
        }
        continue;
      }
      if (isCurrency(v)) continue;               // never route through currency
      const dv = dist.get(v);
      if (dv === undefined || depth + 1 + dv > maxDepth) continue; // cannot reach target in time
      queue.push([...path, v]);
    }
  }

  // Assemble pathway records with per-step reaction alternatives and presence masks.
  const nAccs = meta.accs.length;
  const fullMask = (1n << BigInt(nAccs)) - 1n;
  const pathways = found.map(seq => {
    const steps = [];
    let pathMask = fullMask;
    for (let i = 0; i < seq.length - 1; i++) {
      const rxns = adj.get(seq[i]).get(seq[i + 1]);
      let stepMask = 0n;
      for (const r of rxns) stepMask |= r.mask;
      pathMask &= stepMask;
      steps.push({ from: seq[i], to: seq[i + 1], rxns });
    }
    return { mets: seq, steps, len: seq.length - 1, mask: pathMask, carriers: popcount(pathMask) };
  });
  pathways.sort((a, b) => a.len - b.len || b.carriers - a.carriers);

  // Complete means every simple path up to maxDepth was enumerated; longer
  // pathways are never enumerated and the surface says so either way.
  const complete = !(hitPaths || hitBudget);
  return {
    pathways,
    termination: { complete, hitPaths, hitBudget, expansions, maxDepth, maxPaths, maxExpansions },
  };
}

export function popcount(mask) {
  let c = 0;
  while (mask) { c += Number(mask & 1n); mask >>= 1n; }
  return c;
}

// Per-species carrier breakdown for a pathway mask.
export function carriersBySpecies(mask, meta) {
  const out = new Map();       // species name -> {present: [acc...], total, cells: [{acc, on}]}
  meta.accs.forEach((a, i) => {
    if (!out.has(a.sp)) out.set(a.sp, { present: 0, total: 0, cells: [] });
    const rec = out.get(a.sp);
    const on = ((mask >> BigInt(i)) & 1n) === 1n;
    rec.total++;
    if (on) rec.present++;
    rec.cells.push({ acc: a.acc, on });
  });
  return out;
}
