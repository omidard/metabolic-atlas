// Constraint-based analysis engine: reaction knockouts, growth coupling,
// shadow-price and reduced-cost estimates, production envelopes. Pure compute,
// no DOM. Every number returned here is the result of a real GLPK solve; an
// unsolved or infeasible quantity is returned as null with its solver status,
// never as 0.
//
// The vendored glpk.js 4.0.1 does not expose LP dual values, so shadow prices
// are ESTIMATED by one-sided finite differences: the metabolite's mass-balance
// right-hand side is moved from 0 to +eps and the LP re-solved; the estimate is
// (z(eps) - z(0)) / eps. Within the optimal basis's stability range this equals
// the exact dual; at a degenerate optimum the left and right derivatives can
// differ and the estimate is the right derivative. Reduced costs are derived
// from the estimated duals as rc_j = c_j - sum_i y_i S_ij over the reaction's
// own metabolites. Both tables carry this label in the UI.

import { getGLPK, buildLP, boundType, solveLP, statusName } from './fba.js';

export const ZERO_MU = 1e-6;      // h-1: at or below this a knockout counts as lethal
export const FLUX_TOL = 1e-6;     // mmol gDW-1 h-1: a floor above this counts as coupled
export const COSTLY_FRAC = 0.95;  // mu below this fraction of reference = costly
const BIG = 1e30;

// ---------------------------------------------------------------- target ----
// Product target for coupling-style questions: the guaranteed EXPORT of the
// product. When the GEM carries an exchange for the extracellular form of the
// metabolite's base (whatever compartment form was picked), that exchange is
// the target: its minimum is the guaranteed secretion. Otherwise an added
// demand is the only sink for the metabolite, so its minimum is the guaranteed
// net production. Returns null when the GEM has neither the metabolite nor the
// exchange.
export function analysisTarget(gem, mid) {
  const base = mid.replace(/_[a-z]+$/, '');
  const exId = `EX_${base}_e`;
  if (gem.reactions.some(r => r.ex && r.id === exId)) {
    return { kind: 'exchange', id: exId, extraCols: [], mid };
  }
  if (!gem.metabolites.some(m => m.id === mid)) return null;
  return {
    kind: 'demand', id: `DM_${mid}`, mid,
    extraCols: [{ name: `DM_${mid}`, stoich: { [mid]: -1 }, lb: 0, ub: 1000 }],
  };
}

// --------------------------------------------------------------- session ----
// A session wraps one LP (GEM + medium + optional target column + optional
// permanent knockouts), built once via the shared buildLP so GEM-browser bound
// edits are respected, with index maps for fast bound/row mutation between
// solves. All analysis functions below run on a session.
export async function makeSession(gem, acc, mediumBounds, target, opts = {}) {
  const glpk = await getGLPK();
  const lp = buildLP(glpk, gem, mediumBounds, {
    acc,
    extraCols: target ? target.extraCols : [],
    objective: { direction: 'max', vars: [] },
  });
  const bIdx = new Map(lp.bounds.map(b => [b.name, b]));
  const rowIdx = new Map();
  for (const row of lp.subjectTo) rowIdx.set(row.name, row);
  const S = {
    glpk, gem, acc, mediumBounds, lp, bIdx, rowIdx, target,
    biomass: gem.stats.biomass_id,
    rxnById: new Map(gem.reactions.map(r => [r.id, r])),
    kos: new Set(),
    solves: 0,
  };
  for (const rid of opts.kos || []) applyKO(S, rid);
  return S;
}

export function applyKO(S, rid) {
  const b = S.bIdx.get(rid);
  if (!b || S.kos.has(rid)) return false;
  S.kos.add(rid);
  if (!S._koSaved) S._koSaved = new Map();
  S._koSaved.set(rid, { type: b.type, lb: b.lb, ub: b.ub });
  b.lb = 0; b.ub = 0; b.type = S.glpk.GLP_FX;
  return true;
}

export function removeKO(S, rid) {
  if (!S.kos.has(rid)) return false;
  const saved = S._koSaved.get(rid);
  Object.assign(S.bIdx.get(rid), saved);
  S._koSaved.delete(rid);
  S.kos.delete(rid);
  return true;
}

function saveBound(S, name) {
  const b = S.bIdx.get(name);
  return b ? { name, type: b.type, lb: b.lb, ub: b.ub } : null;
}
function restoreBound(S, saved) {
  if (saved) Object.assign(S.bIdx.get(saved.name), { type: saved.type, lb: saved.lb, ub: saved.ub });
}
function setFixed(S, name, val) {
  const b = S.bIdx.get(name);
  b.lb = val; b.ub = val; b.type = S.glpk.GLP_FX;
}
function setRange(S, name, lb, ub) {
  const b = S.bIdx.get(name);
  b.lb = lb; b.ub = ub; b.type = boundType(S.glpk, lb, ub);
}

async function opt(S, dir, varName) {
  S.lp.objective = {
    direction: dir === 'min' ? S.glpk.GLP_MIN : S.glpk.GLP_MAX,
    name: 'obj',
    vars: [{ name: varName, coef: 1 }],
  };
  S.solves++;
  return solveLP(S.glpk, S.lp);
}

// ------------------------------------------------------------ base state ----
// Reference state of the session's model: max growth, and when a target is
// set, max product (growth free) and the guaranteed product floor at biomass
// held at >= frac of max growth. The floor is the coupling test: floor >
// FLUX_TOL means growth-coupled at that fraction.
export async function baseState(S, frac) {
  const out = { frac, solves0: S.solves };
  const g = await opt(S, 'max', S.biomass);
  out.mu = g.optimal ? g.z : null;
  out.muStatus = g.status;
  out.muVars = g.optimal ? g.vars : null;
  if (S.target) {
    const p = await opt(S, 'max', S.target.id);
    out.productMax = p.optimal ? p.z : null;
    out.productStatus = p.status;
    if (out.mu != null && out.mu > ZERO_MU) {
      const saved = saveBound(S, S.biomass);
      setRange(S, S.biomass, frac * out.mu, Math.max(saved.ub, out.mu));
      const f = await opt(S, 'min', S.target.id);
      restoreBound(S, saved);
      out.floor = f.optimal ? f.z : null;
      out.floorStatus = f.status;
      out.floorVars = f.optimal ? f.vars : null;
      out.coupled = out.floor != null ? out.floor > FLUX_TOL : null;
    } else {
      out.floor = null;
      out.coupled = null;
    }
  }
  out.solves = S.solves - out.solves0;
  return out;
}

// -------------------------------------------------------------- KO sweep ----
// Single-reaction knockouts over every non-exchange, non-biomass reaction.
// Per knockout: max biomass; when a target is set and the knockout is not
// lethal, also max product (growth free) and the product floor at biomass >=
// frac * the knockout's own max. A reaction whose current bounds are already
// [0, 0] is reported from the reference solve (the LP is identical; no solve
// is spent) and marked so. cbRow fires after every finished row.
export function sweepScope(gem, biomassId) {
  const rxns = gem.reactions.filter(r => !r.ex && r.id !== biomassId);
  return { rxns, excluded: gem.reactions.length - rxns.length, total: gem.reactions.length };
}

export async function koSweep(S, refs, opts = {}) {
  const { rxns, excluded, total } = sweepScope(S.gem, S.biomass);
  const frac = refs.frac;
  const rows = [];
  let cancelled = false;
  const solves0 = S.solves;
  for (let i = 0; i < rxns.length; i++) {
    if (opts.shouldStop && opts.shouldStop()) { cancelled = true; break; }
    const r = rxns[i];
    const b = S.bIdx.get(r.id);
    let row;
    if (b.lb === 0 && b.ub === 0) {
      row = {
        id: r.id, name: r.name || '', genes: r.genes || [],
        alreadyClosed: true, lethalNoSteadyState: false,
        mu: refs.mu, productMax: refs.productMax ?? null, floor: refs.floor ?? null,
      };
    } else {
      const saved = saveBound(S, r.id);
      setFixed(S, r.id, 0);
      const g = await opt(S, 'max', S.biomass);
      row = {
        id: r.id, name: r.name || '', genes: r.genes || [],
        alreadyClosed: false,
        mu: g.optimal ? g.z : null,
        muStatus: g.status,
        lethalNoSteadyState: !g.optimal,
        productMax: null, floor: null,
      };
      const muKO = row.mu;
      if (S.target && muKO != null && muKO > ZERO_MU) {
        const p = await opt(S, 'max', S.target.id);
        row.productMax = p.optimal ? p.z : null;
        const savedBio = saveBound(S, S.biomass);
        setRange(S, S.biomass, frac * muKO, Math.max(savedBio.ub, muKO));
        const f = await opt(S, 'min', S.target.id);
        restoreBound(S, savedBio);
        row.floor = f.optimal ? f.z : null;
      }
      restoreBound(S, saved);
    }
    row.cls = classifyKO(row, refs);
    rows.push(row);
    if (opts.cbRow) opts.cbRow(row, i + 1, rxns.length);
  }
  return {
    rows, swept: rows.length, scope: rxns.length, excluded, total,
    cancelled, solves: S.solves - solves0,
  };
}

// Classification thresholds are stated on the methods page. A knockout can be
// both costly and floor-raising; the floor takes precedence in the label and
// both numbers stay visible in the table.
export function classifyKO(row, refs) {
  const muRef = refs.mu;
  if (row.lethalNoSteadyState || (row.mu != null && row.mu <= ZERO_MU)) return 'lethal';
  if (refs.floor != null && row.floor != null && row.floor > Math.max(refs.floor, 0) + FLUX_TOL) return 'beneficial';
  if (muRef != null && row.mu != null && row.mu < COSTLY_FRAC * muRef) return 'costly';
  if (row.mu == null) return 'unsolved';
  return 'neutral';
}

// -------------------------------------------- shadow prices, reduced costs ----
// Finite-difference shadow-price estimates at the max-growth optimum (see the
// module header for the method and its degeneracy caveat). eps defaults to
// 1e-3 mmol gDW-1 h-1. Results are cached on the session per metabolite.
export const SHADOW_EPS = 1e-3;

export async function shadowPrices(S, metIds, opts = {}) {
  const eps = opts.eps || SHADOW_EPS;
  if (!S._dualCache) S._dualCache = new Map();
  if (S._dualZ0 == null) {
    const g = await opt(S, 'max', S.biomass);
    if (!g.optimal) return { ok: false, status: g.status, statusText: statusName(g.status) };
    S._dualZ0 = g.z;
  }
  const out = new Map();
  let done = 0;
  for (const mid of metIds) {
    if (opts.shouldStop && opts.shouldStop()) return { ok: true, z0: S._dualZ0, values: out, cancelled: true };
    if (S._dualCache.has(mid)) { out.set(mid, S._dualCache.get(mid)); done++; continue; }
    const row = S.rowIdx.get(mid);
    let rec;
    if (!row) {
      rec = { y: null, note: 'not in this GEM' };
    } else {
      const saved = { ...row.bnds };
      row.bnds = { type: S.glpk.GLP_FX, lb: eps, ub: eps };
      let s = await opt(S, 'max', S.biomass);
      if (s.optimal) {
        rec = { y: (s.z - S._dualZ0) / eps, side: 'right' };
      } else {
        // +eps infeasible (e.g. a conserved cofactor pool cannot carry a
        // surplus): try the left derivative at -eps before giving up.
        row.bnds = { type: S.glpk.GLP_FX, lb: -eps, ub: -eps };
        s = await opt(S, 'max', S.biomass);
        rec = s.optimal
          ? { y: (S._dualZ0 - s.z) / eps, side: 'left' }
          : { y: null, note: 'not finite by perturbation (both signs infeasible; conserved pool or blocked metabolite)' };
      }
      row.bnds = saved;
    }
    S._dualCache.set(mid, rec);
    out.set(mid, rec);
    done++;
    if (opts.onProgress) opts.onProgress(done, metIds.length);
  }
  return { ok: true, z0: S._dualZ0, values: out, cancelled: false };
}

// Reduced-cost estimates for a set of reactions, derived from shadow-price
// estimates of exactly the metabolites those reactions touch (computed or
// reused from cache). c_j is 1 for the biomass reaction and 0 otherwise.
export async function reducedCosts(S, rxnIds, opts = {}) {
  const mets = new Set();
  const rxns = [];
  for (const rid of rxnIds) {
    const r = S.rxnById.get(rid);
    if (!r) continue;
    rxns.push(r);
    for (const m of Object.keys(r.stoich)) mets.add(m);
  }
  const sp = await shadowPrices(S, [...mets], opts);
  if (!sp.ok) return sp;
  const out = new Map();
  for (const r of rxns) {
    let rc = r.id === S.biomass ? 1 : 0;
    let missing = null;
    for (const [m, c] of Object.entries(r.stoich)) {
      const rec = sp.values.get(m);
      if (!rec || rec.y == null) { missing = m; break; }
      rc -= rec.y * c;
    }
    out.set(r.id, missing ? { rc: null, note: `needs shadow price of ${missing}` } : { rc });
  }
  return { ok: true, z0: sp.z0, values: out, cancelled: sp.cancelled };
}

// ------------------------------------------------------- coupling search ----
// Greedy knockout search for growth coupling. Definitions:
//   mu    = max biomass with the current knockout set
//   floor = min product at biomass >= frac * mu   (coupled when floor > tol)
//   mu0   = max biomass with product export capped at the coupling tolerance
//           (the decoupled optimum; a cap, not an exact zero fix, so the
//           decoupled solve and the floor test share one tolerance)
// Coupling at fraction frac is equivalent to mu0 < frac * mu, so the greedy
// score of a candidate knockout is the ratio mu0/mu it produces; the search
// applies the candidate with the lowest ratio and repeats. Candidates are the
// reactions carrying flux in a parsimonious zero-product max-growth solution:
// a knockout that leaves that solution intact cannot lower mu0, so only its
// support can help. This is a heuristic: it can miss knockout sets a global
// search would find, and a failure to couple within K knockouts is reported as
// exactly that, never as impossibility.
export async function couplingSearch(S, opts = {}) {
  const frac = opts.frac ?? 0.9;
  const K = opts.K ?? 3;
  const requireGene = opts.requireGene !== false;
  const log = opts.onIteration || (() => {});
  const solves0 = S.solves;
  const iterations = [];
  const koSet = [];

  const fin = (verdict, reason, extra = {}) => ({
    verdict, reason, koSet, iterations, frac, K,
    solves: S.solves - solves0, ...extra,
  });

  // base viability and capability
  let g = await opt(S, 'max', S.biomass);
  if (!g.optimal || g.z <= ZERO_MU) {
    return fin('no-growth',
      `The model's max growth on this medium is ${g.optimal ? g.z.toExponential(2) : statusName(g.status)}; a coupling search needs a growing reference state.`);
  }
  let mu = g.z;
  const minMu = opts.minMu ?? 0.1 * mu;
  const pmaxS = await opt(S, 'max', S.target.id);
  const productMax = pmaxS.optimal ? pmaxS.z : null;
  if (productMax == null || productMax <= FLUX_TOL) {
    return fin('no-production',
      'The product cannot be exported at all on this GEM and medium (max product flux at or below tolerance). Knockouts only remove capability, so no knockout set can create the missing route.', { mu, productMax });
  }

  for (let it = 0; it <= K; it++) {
    if (opts.shouldStop && opts.shouldStop()) return fin('cancelled', 'Search cancelled.', { mu });

    // coupling test at the current knockout set
    const savedBio = saveBound(S, S.biomass);
    setRange(S, S.biomass, frac * mu, Math.max(savedBio.ub, mu));
    const f = await opt(S, 'min', S.target.id);
    restoreBound(S, savedBio);
    const floor = f.optimal ? f.z : null;
    if (floor != null && floor > FLUX_TOL) {
      return fin('coupled', `Guaranteed product flux at biomass >= ${(frac * 100).toFixed(0)}% of max growth is positive.`,
        { mu, floor, floorVars: f.vars, productMax });
    }

    // decoupled optimum mu0 and its ratio
    const savedT = saveBound(S, S.target.id);
    setRange(S, S.target.id, Math.min(savedT.lb, 0), FLUX_TOL);
    const g0 = await opt(S, 'max', S.biomass);
    const mu0 = g0.optimal ? g0.z : 0;   // infeasible: growth without product impossible
    const ratio = mu0 / mu;

    if (it === K) {
      restoreBound(S, savedT);
      return fin('not-coupled',
        `${K} knockout${K === 1 ? '' : 's'} did not couple the product at the ${(frac * 100).toFixed(0)}% growth fraction (decoupled-growth ratio ${ratio.toFixed(3)}; coupling needs it below ${frac}). A larger knockout set or a different medium may; this greedy single-knockout search does not prove impossibility.`,
        { mu, floor, mu0, ratio, productMax });
    }

    // candidates: support of a parsimonious zero-product max-growth solution
    const support = await pfbaSupport(S, mu0);
    restoreBound(S, savedT);
    const candidates = [];
    for (const rid of support) {
      const r = S.rxnById.get(rid);
      if (!r || r.ex || rid === S.biomass || S.kos.has(rid)) continue;
      if (rid.startsWith('DM_') || rid.startsWith('SK_')) continue;
      if (requireGene && !(r.gpr && r.gpr.trim())) continue;
      candidates.push(rid);
    }
    if (!candidates.length) {
      return fin('not-coupled',
        'No eligible candidate knockout carries flux in the zero-product optimum' + (requireGene ? ' (reactions without a gene rule are excluded)' : '') + ', so no single knockout can lower decoupled growth.',
        { mu, floor, mu0, ratio, productMax });
    }

    // score each candidate by the mu0/mu ratio it produces
    let best = null;
    let tested = 0, skippedLethal = 0;
    for (const rid of candidates) {
      if (opts.shouldStop && opts.shouldStop()) return fin('cancelled', 'Search cancelled.', { mu });
      const saved = saveBound(S, rid);
      setFixed(S, rid, 0);
      const gc = await opt(S, 'max', S.biomass);
      const muC = gc.optimal ? gc.z : null;
      if (muC == null || muC < minMu) {
        skippedLethal++;
        restoreBound(S, saved);
        continue;
      }
      const sT = saveBound(S, S.target.id);
      setRange(S, S.target.id, Math.min(sT.lb, 0), FLUX_TOL);
      const g0c = await opt(S, 'max', S.biomass);
      restoreBound(S, sT);
      restoreBound(S, saved);
      const mu0C = g0c.optimal ? g0c.z : 0;
      const ratioC = mu0C / muC;
      tested++;
      if (!best || ratioC < best.ratio - 1e-12) best = { rid, mu: muC, mu0: mu0C, ratio: ratioC };
      if (opts.onCandidate) opts.onCandidate(tested, candidates.length, it + 1);
    }
    iterations.push({
      iteration: it + 1, mu, mu0, ratio, floor,
      candidates: candidates.length, tested, skippedLethal,
      best: best ? { id: best.rid, ratio: best.ratio } : null,
    });
    log(iterations[iterations.length - 1]);

    if (!best || best.ratio >= ratio - 1e-9) {
      return fin('not-coupled',
        `None of the ${tested} viable candidate knockouts lowered the decoupled-growth ratio below the current ${ratio.toFixed(3)} (${skippedLethal} of ${candidates.length} candidates dropped max growth under the viability floor).`,
        { mu, floor, mu0, ratio, productMax });
    }

    applyKO(S, best.rid);
    const rec = S.rxnById.get(best.rid);
    koSet.push({ id: best.rid, name: rec.name || '', genes: rec.genes || [], gpr: rec.gpr || '', ratioAfter: best.ratio, muAfter: best.mu });
    mu = best.mu;
  }
  return fin('not-coupled', 'Search ended without coupling.', { mu });
}

// Support (|v| > 1e-6) of a parsimonious solution at biomass >= 0.999 * muFix
// under the session's CURRENT bounds (callers fix the product beforehand).
// Built as a fresh LP with |v| helper columns, like the Mode-2 pFBA.
async function pfbaSupport(S, muFix) {
  if (muFix <= ZERO_MU) return [];
  const glpk = S.glpk;
  const lp = {
    name: 'pfba_support',
    objective: { direction: glpk.GLP_MIN, name: 'total_flux', vars: [] },
    subjectTo: S.lp.subjectTo.map(r => ({ name: r.name, vars: r.vars, bnds: r.bnds })),
    bounds: S.lp.bounds.map(b => ({ ...b })),
  };
  for (const b of lp.bounds) if (b.name === S.biomass) {
    b.lb = 0.999 * muFix;
    b.ub = Math.max(b.ub, muFix);
    b.type = boundType(glpk, b.lb, b.ub);
  }
  const absVars = [];
  const extra = [];
  for (const r of S.gem.reactions) {
    const a = 'abs_' + r.id;
    absVars.push({ name: a, coef: 1 });
    lp.bounds.push({ name: a, type: glpk.GLP_LO, lb: 0, ub: BIG });
    extra.push({ name: 'ap_' + r.id, vars: [{ name: a, coef: 1 }, { name: r.id, coef: -1 }], bnds: { type: glpk.GLP_LO, lb: 0, ub: 0 } });
    extra.push({ name: 'an_' + r.id, vars: [{ name: a, coef: 1 }, { name: r.id, coef: 1 }], bnds: { type: glpk.GLP_LO, lb: 0, ub: 0 } });
  }
  lp.subjectTo = lp.subjectTo.concat(extra);
  lp.objective.vars = absVars;
  S.solves++;
  const s = await solveLP(glpk, lp);
  if (!s.optimal) {
    // fall back to the support of a plain (non-parsimonious) solve
    const g = await opt(S, 'max', S.biomass);
    if (!g.optimal) return [];
    return S.gem.reactions.filter(r => Math.abs(g.vars[r.id] ?? 0) > 1e-6).map(r => r.id);
  }
  return S.gem.reactions.filter(r => Math.abs(s.vars[r.id] ?? 0) > 1e-6).map(r => r.id);
}

// ---------------------------------------------------- production envelope ----
// Max and min product flux at each of nPoints biomass levels from 0 to max
// growth, under the session's current bounds and knockouts. Points where the
// LP fails are returned with null bounds and the solver status.
export async function productionEnvelope(S, nPoints = 24, opts = {}) {
  const solves0 = S.solves;
  const g = await opt(S, 'max', S.biomass);
  if (!g.optimal) return { ok: false, status: g.status, statusText: statusName(g.status) };
  const mu = g.z;
  const points = [];
  const saved = saveBound(S, S.biomass);
  for (let k = 0; k < nPoints; k++) {
    if (opts.shouldStop && opts.shouldStop()) break;
    let b = mu * k / (nPoints - 1);
    if (k === nPoints - 1) b = mu;
    setFixed(S, S.biomass, b);
    let mx = await opt(S, 'max', S.target.id);
    let mn = mx.optimal ? await opt(S, 'min', S.target.id) : mx;
    if (!mx.optimal && k === nPoints - 1) {
      // numerical top point: back off by 1e-6 relative
      setFixed(S, S.biomass, b * (1 - 1e-6));
      mx = await opt(S, 'max', S.target.id);
      mn = mx.optimal ? await opt(S, 'min', S.target.id) : mx;
    }
    points.push({
      mu: b,
      max: mx.optimal ? mx.z : null,
      min: mn.optimal ? mn.z : null,
      status: mx.optimal && mn.optimal ? null : statusName((mx.optimal ? mn : mx).status),
    });
    if (opts.onProgress) opts.onProgress(k + 1, nPoints);
  }
  restoreBound(S, saved);
  return { ok: true, mu, points, requested: nPoints, solves: S.solves - solves0 };
}

// ------------------------------------------------------------------ yield ----
// Yield of the product floor against substrate uptake in the floor solution.
// mmol/mmol is the flux ratio; C-mol/C-mol additionally needs both formulas.
export function carbonCount(formula) {
  if (!formula) return null;
  const m = /(?:^|[^A-Za-z])C([0-9]*)(?![a-z])/.exec(formula);
  if (!m) return null;
  return m[1] ? parseInt(m[1], 10) : 1;
}

export function floorYield(S, floorVars, floor, substrateEx, substrateMid, productMid) {
  if (!floorVars || floor == null || !substrateEx) {
    return { mmol: null, cmol: null, note: 'not computed (no substrate uptake identified)' };
  }
  const up = floorVars[substrateEx];
  if (up == null || up >= -FLUX_TOL) {
    return { mmol: null, cmol: null, note: 'not computed (no substrate uptake in the floor solution)' };
  }
  const mmol = floor / Math.abs(up);
  const fSub = S.gem.metabolites.find(m => m.id === substrateMid);
  const fProd = S.gem.metabolites.find(m => m.id === productMid);
  const cS = carbonCount(fSub && fSub.formula);
  const cP = carbonCount(fProd && fProd.formula);
  const cmol = (cS && cP) ? mmol * cP / cS : null;
  return {
    mmol, cmol, uptake: Math.abs(up),
    note: cmol == null ? 'C-mol basis not computed (formula missing a carbon count)' : null,
  };
}
