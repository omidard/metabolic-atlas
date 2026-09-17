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

// Parsimonious flux distribution (min sum |v|) under the session's CURRENT
// bounds, with optional per-column bound overrides applied to a fresh LP copy;
// the session LP itself is untouched. Built with |v| helper columns, like the
// Mode-2 pFBA. Returns the flux of every session column (model reactions plus
// any target column), or {optimal:false} with the solver status.
export async function pfbaFluxes(S, opts = {}) {
  const glpk = S.glpk;
  const lp = {
    name: 'pfba',
    objective: { direction: glpk.GLP_MIN, name: 'total_flux', vars: [] },
    subjectTo: S.lp.subjectTo.map(r => ({ name: r.name, vars: r.vars, bnds: r.bnds })),
    bounds: S.lp.bounds.map(b => ({ ...b })),
  };
  for (const f of opts.fix || []) {
    for (const b of lp.bounds) if (b.name === f.name) {
      b.lb = f.lb; b.ub = f.ub;
      b.type = boundType(glpk, f.lb, f.ub);
    }
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
  if (!s.optimal) return { optimal: false, status: s.status };
  const fluxes = {};
  for (const b of S.lp.bounds) fluxes[b.name] = s.vars[b.name] ?? null;
  return { optimal: true, fluxes, totalFlux: s.z };
}

// Support (|v| > 1e-6) of a parsimonious solution at biomass >= 0.999 * muFix
// under the session's CURRENT bounds (callers fix the product beforehand).
async function pfbaSupport(S, muFix) {
  if (muFix <= ZERO_MU) return [];
  const bio = S.bIdx.get(S.biomass);
  const p = await pfbaFluxes(S, {
    fix: [{ name: S.biomass, lb: 0.999 * muFix, ub: Math.max(bio.ub, muFix) }],
  });
  if (!p.optimal) {
    // fall back to the support of a plain (non-parsimonious) solve
    const g = await opt(S, 'max', S.biomass);
    if (!g.optimal) return [];
    return S.gem.reactions.filter(r => Math.abs(g.vars[r.id] ?? 0) > 1e-6).map(r => r.id);
  }
  return S.gem.reactions.filter(r => Math.abs(p.fluxes[r.id] ?? 0) > 1e-6).map(r => r.id);
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

// ------------------------------------------------------------------ FSEOF ----
// Flux scanning with enforced objective flux: the product export is enforced
// at nSteps+1 evenly spaced levels from 0 to maxFrac * max product; at each
// level biomass is maximised and the flux distribution is made parsimonious
// (pFBA at >= 99.9% of that level's own max growth). A reaction whose |flux|
// rises monotonically (within tol) across the solved levels by more than
// FSEOF_MIN_CHANGE is an amplification (over-expression) target; one whose
// |flux| falls monotonically is an attenuation (down-regulation) target.
// Reactions that change flux sign across the scan, or are never active, are
// classified as neither. The slope is the least-squares slope of |flux|
// against the enforced product flux. Exchanges, the biomass reaction and the
// target column are out of scope. The top enforced level defaults to 95% of
// the max product because at 100% the optimum is a single point where the
// biomass LP is often degenerate.
export const FSEOF_STEPS = 10;
export const FSEOF_MAX_FRAC = 0.95;
export const FSEOF_MIN_CHANGE = 1e-4;   // mmol gDW-1 h-1 of |flux| change across the scan
const FSEOF_TOL = 1e-6;                 // monotonicity slack per step

export async function fseofScan(S, opts = {}) {
  const nSteps = opts.nSteps ?? FSEOF_STEPS;
  const maxFrac = opts.maxFrac ?? FSEOF_MAX_FRAC;
  const solves0 = S.solves;
  const fin = (o) => ({ nSteps, maxFrac, solves: S.solves - solves0, ...o });
  if (!S.target) return fin({ ok: false, reason: 'no product target is set' });

  const g = await opt(S, 'max', S.biomass);
  if (!g.optimal || g.z <= ZERO_MU) {
    return fin({
      ok: false,
      reason: `max growth on this medium is ${g.optimal ? g.z.toExponential(2) + ' 1/h' : statusName(g.status)}; the scan maximises biomass at each enforced level and needs a growing base state`,
    });
  }
  const p = await opt(S, 'max', S.target.id);
  if (!p.optimal || p.z <= FLUX_TOL) {
    return fin({
      ok: false, productMax: p.optimal ? p.z : null,
      reason: 'the product cannot be exported on this GEM and medium (max product flux at or below tolerance), so there is no objective flux to enforce',
    });
  }
  const pmax = p.z;
  const fMax = maxFrac * pmax;
  const scope = S.gem.reactions.filter(r => !r.ex && r.id !== S.biomass && r.id !== S.target.id);
  const saved = saveBound(S, S.target.id);
  const bio = S.bIdx.get(S.biomass);
  const levels = [];
  const series = new Map(scope.map(r => [r.id, []]));
  let cancelled = false;
  for (let k = 0; k <= nSteps; k++) {
    if (opts.shouldStop && opts.shouldStop()) { cancelled = true; break; }
    const f = fMax * k / nSteps;
    setRange(S, S.target.id, f, Math.max(saved.ub, fMax));
    const gk = await opt(S, 'max', S.biomass);
    if (!gk.optimal) {
      levels.push({ f, mu: null, solved: false, status: gk.status });
      continue;
    }
    const pk = await pfbaFluxes(S, {
      fix: [{ name: S.biomass, lb: 0.999 * gk.z, ub: Math.max(bio.ub, gk.z) }],
    });
    const fl = pk.optimal ? pk.fluxes : gk.vars;
    levels.push({ f, mu: gk.z, solved: true, pfba: pk.optimal });
    for (const r of scope) series.get(r.id).push({ f, v: fl[r.id] ?? null });
    if (opts.onProgress) opts.onProgress(k + 1, nSteps + 1);
  }
  restoreBound(S, saved);

  const solvedLevels = levels.filter(l => l.solved);
  if (solvedLevels.length < 3) {
    return fin({
      ok: false, productMax: pmax, fMax, levels, cancelled,
      reason: `only ${solvedLevels.length} of ${levels.length} enforced levels solved; at least 3 are needed to read a trend`,
    });
  }

  const up = [], down = [];
  let active = 0, signChanging = 0;
  for (const r of scope) {
    const pts = series.get(r.id).filter(q => q.v != null);
    if (pts.length !== solvedLevels.length) continue;   // a level lacked this flux; unclassifiable
    const a = pts.map(q => Math.abs(q.v));
    if (Math.max(...a) < FSEOF_TOL) continue;           // never active in the scan
    active++;
    const pos = pts.some(q => q.v > FSEOF_TOL);
    const neg = pts.some(q => q.v < -FSEOF_TOL);
    if (pos && neg) { signChanging++; continue; }       // direction flips; neither list
    let inc = true, dec = true;
    for (let i = 1; i < a.length; i++) {
      if (a[i] < a[i - 1] - FSEOF_TOL) inc = false;
      if (a[i] > a[i - 1] + FSEOF_TOL) dec = false;
    }
    const rise = a[a.length - 1] - a[0];
    // least-squares slope of |v| against enforced product flux
    const n = pts.length;
    const mx = pts.reduce((s2, q) => s2 + q.f, 0) / n;
    const my = a.reduce((s2, v) => s2 + v, 0) / n;
    let sxy = 0, sxx = 0;
    for (let i = 0; i < n; i++) { sxy += (pts[i].f - mx) * (a[i] - my); sxx += (pts[i].f - mx) ** 2; }
    const slope = sxx > 0 ? sxy / sxx : 0;
    const rec = {
      id: r.id, name: r.name || '', genes: r.genes || [], gpr: r.gpr || '',
      subsystem: r.subsystem || '', slope,
      v0: pts[0].v, vEnd: pts[n - 1].v, absChange: rise,
      fluxes: pts.map(q => q.v),
    };
    if (inc && rise > FSEOF_MIN_CHANGE) up.push(rec);
    else if (dec && -rise > FSEOF_MIN_CHANGE) down.push(rec);
  }
  up.sort((x, y) => y.slope - x.slope);
  down.sort((x, y) => x.slope - y.slope);
  return fin({
    ok: true, mu: g.z, productMax: pmax, fMax, levels, cancelled,
    up, down, scanned: scope.length, activeInScan: active, signChanging,
    excluded: S.gem.reactions.length - scope.length,
    total: S.gem.reactions.length,
    solvedLevels: solvedLevels.length,
  });
}

// ----------------------------------------------------------- linear MOMA ----
// Wild-type reference for MOMA: the parsimonious flux distribution at >= 99.9%
// of max growth under the session's current bounds. pFBA picks ONE of the
// possibly many optimal distributions; the MOMA prediction is relative to it.
export async function wtReference(S) {
  const g = await opt(S, 'max', S.biomass);
  if (!g.optimal || g.z <= ZERO_MU) {
    return { ok: false, mu: g.optimal ? g.z : null, status: g.status, statusText: statusName(g.status) };
  }
  const bio = S.bIdx.get(S.biomass);
  const p = await pfbaFluxes(S, {
    fix: [{ name: S.biomass, lb: 0.999 * g.z, ub: Math.max(bio.ub, g.z) }],
  });
  if (!p.optimal) return { ok: false, mu: g.z, status: p.status, statusText: statusName(p.status) };
  return { ok: true, mu: g.z, fluxes: p.fluxes, totalFlux: p.totalFlux };
}

// Linear (L1) MOMA: minimise sum |v - v_wt| subject to S.v = 0 and the session
// bounds with the given reactions knocked out. Each deviation splits into
// dp - dn with dp, dn >= 0 and the row v - dp + dn = v_wt, so the LP minimises
// the exact L1 distance. This is the linear variant; the quadratic (L2) MOMA
// objective is not solvable with the LP-only solver in this build.
export async function linearMOMA(S, wtFluxes, koIds, opts = {}) {
  const glpk = S.glpk;
  const solves0 = S.solves;
  const lp = {
    name: 'lin_moma',
    objective: { direction: glpk.GLP_MIN, name: 'l1_dist', vars: [] },
    subjectTo: S.lp.subjectTo.map(r => ({ name: r.name, vars: r.vars, bnds: r.bnds })),
    bounds: S.lp.bounds.map(b => ({ ...b })),
  };
  const missing = [];
  for (const rid of koIds) {
    const b = lp.bounds.find(x => x.name === rid);
    if (!b) { missing.push(rid); continue; }
    b.lb = 0; b.ub = 0; b.type = glpk.GLP_FX;
  }
  const devVars = [];
  const extra = [];
  for (const b of S.lp.bounds) {
    const wt = wtFluxes[b.name];
    if (wt == null) continue;
    devVars.push({ name: 'dp_' + b.name, coef: 1 }, { name: 'dn_' + b.name, coef: 1 });
    lp.bounds.push(
      { name: 'dp_' + b.name, type: glpk.GLP_LO, lb: 0, ub: BIG },
      { name: 'dn_' + b.name, type: glpk.GLP_LO, lb: 0, ub: BIG });
    extra.push({
      name: 'dev_' + b.name,
      vars: [{ name: b.name, coef: 1 }, { name: 'dp_' + b.name, coef: -1 }, { name: 'dn_' + b.name, coef: 1 }],
      bnds: { type: glpk.GLP_FX, lb: wt, ub: wt },
    });
  }
  lp.subjectTo = lp.subjectTo.concat(extra);
  lp.objective.vars = devVars;
  S.solves++;
  const s = await solveLP(glpk, lp);
  if (!s.optimal) {
    return { ok: false, status: s.status, statusText: statusName(s.status), missing, solves: S.solves - solves0 };
  }
  const fluxes = {};
  for (const b of S.lp.bounds) fluxes[b.name] = s.vars[b.name] ?? null;
  return {
    ok: true, distance: s.z, fluxes,
    mu: fluxes[S.biomass],
    product: S.target ? fluxes[S.target.id] : null,
    missing, solves: S.solves - solves0,
  };
}

// FBA prediction for the same knockout set, for the side-by-side comparison:
// max growth with the knockouts applied, and the product flux in a
// parsimonious solution at that optimum (one of possibly many optima).
export async function fbaKnockout(S, koIds) {
  const saved = [];
  for (const rid of koIds) {
    const sb = saveBound(S, rid);
    if (sb) { saved.push(sb); setFixed(S, rid, 0); }
  }
  const g = await opt(S, 'max', S.biomass);
  const out = {
    mu: g.optimal ? g.z : null, status: g.status, statusText: statusName(g.status),
    product: null, pfbaOptimal: false,
  };
  if (g.optimal && g.z > ZERO_MU) {
    const bio = S.bIdx.get(S.biomass);
    const p = await pfbaFluxes(S, {
      fix: [{ name: S.biomass, lb: 0.999 * g.z, ub: Math.max(bio.ub, g.z) }],
    });
    if (p.optimal) {
      out.pfbaOptimal = true;
      out.product = S.target ? (p.fluxes[S.target.id] ?? null) : null;
      out.fluxes = p.fluxes;
    }
  }
  for (const sb of saved) restoreBound(S, sb);
  return out;
}

// ---------------------------------------------------------- flux sampling ----
// Random-objective vertex sampling of the feasible flux space: each sample is
// the optimal solution of one random dense linear objective (coefficients
// drawn standard-normal from a seeded generator) over the session's polytope,
// optionally with biomass held at >= biomassFrac of max growth. Every sample
// is therefore a VERTEX of the polytope; the empirical distribution is over
// vertices weighted by the random-direction measure, NOT a uniform sample of
// the feasible space. Objectives whose LP fails the optimality or bound check
// are counted as failed and contribute no sample.
export function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export const SAMPLE_DEFAULT_N = 200;
export const SAMPLE_MAX_N = 500;      // browser cap; stated in the UI

export async function sampleFluxSpace(S, opts = {}) {
  const n = Math.min(opts.n ?? SAMPLE_DEFAULT_N, SAMPLE_MAX_N);
  const frac = opts.biomassFrac ?? 0;
  const seed = opts.seed ?? 1;
  const solves0 = S.solves;
  let mu = null, savedBio = null;
  if (frac > 0) {
    const g = await opt(S, 'max', S.biomass);
    if (!g.optimal || g.z <= ZERO_MU) {
      return {
        ok: false, solves: S.solves - solves0,
        reason: `max growth on this medium is ${g.optimal ? g.z.toExponential(2) + ' 1/h' : statusName(g.status)}; a biomass-fraction constraint needs a growing state (set the fraction to 0 to sample without it)`,
      };
    }
    mu = g.z;
    savedBio = saveBound(S, S.biomass);
    setRange(S, S.biomass, frac * mu, Math.max(savedBio.ub, mu));
  }
  const rids = [...S.gem.reactions.map(r => r.id), ...(S.target && S.target.kind === 'demand' ? [S.target.id] : [])];
  const store = new Map(rids.map(rid => [rid, new Float64Array(n)]));
  const rand = mulberry32(seed);
  // Box-Muller standard normals
  const randn = () => {
    let u = 0, v = 0;
    while (u === 0) u = rand();
    while (v === 0) v = rand();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  };
  let done = 0, failed = 0;
  let cancelled = false;
  for (let i = 0; i < n; i++) {
    if (opts.shouldStop && opts.shouldStop()) { cancelled = true; break; }
    S.lp.objective = {
      direction: S.glpk.GLP_MAX, name: 'rand',
      vars: rids.map(rid => ({ name: rid, coef: randn() })),
    };
    S.solves++;
    const s = await solveLP(S.glpk, S.lp);
    if (!s.optimal) { failed++; continue; }
    for (const rid of rids) store.get(rid)[done] = s.vars[rid] ?? 0;
    done++;
    if (opts.onProgress) opts.onProgress(i + 1, n, done, failed);
  }
  if (savedBio) restoreBound(S, savedBio);

  const quant = (sorted, q) => {
    const t = q * (sorted.length - 1);
    const lo = Math.floor(t), hi = Math.ceil(t);
    return sorted[lo] + (sorted[hi] - sorted[lo]) * (t - lo);
  };
  const stats = new Map();
  if (done > 0) {
    for (const rid of rids) {
      const vals = store.get(rid).slice(0, done);
      const sorted = Float64Array.from(vals).sort();
      const mean = vals.reduce((a, b) => a + b, 0) / done;
      stats.set(rid, {
        mean,
        median: quant(sorted, 0.5),
        p5: quant(sorted, 0.05),
        p95: quant(sorted, 0.95),
        min: sorted[0],
        max: sorted[sorted.length - 1],
      });
    }
  }
  return {
    ok: done > 0, requested: n, samples: done, failed, cancelled,
    mu, frac, seed, stats,
    raw: store, rids,
    sampler: 'random-objective vertex sampling',
    solves: S.solves - solves0,
    reason: done > 0 ? null : `none of the ${n - (cancelled ? n - done - failed : 0)} attempted objectives returned an optimal, bound-respecting solution`,
  };
}

// ------------------------------------------------------------------ yield ----
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
