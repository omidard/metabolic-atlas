// In-browser flux engine for Mode 2 (GLPK 5.0 via the vendored glpk.js 4.0.1
// WASM worker). Builds S.v = 0 with bounds from a GEM json plus a medium, and
// answers: growth on the medium, pathway flux feasibility, pFBA distribution,
// FVA ranges. Same LP pattern as the EcopanGEM Flux Analysis Studio engine.
//
// Medium convention: an exchange listed in mediumBounds gets that lower bound;
// every exchange NOT listed is closed (lb = 0). Upper bounds (secretion) keep
// the model's own values. In-memory bound edits from the GEM browser override
// the model's bounds for non-exchange reactions and the ub of exchanges.

import GLPK from '../vendor/glpk.esm.js';
import { getEdit } from './data.js';

let _glpk = null;
export async function getGLPK() {
  if (!_glpk) _glpk = await GLPK();
  return _glpk;
}

const BIG = 1e30;
export const STEP_MIN_FLUX = 1e-4;   // mmol gDW-1 h-1 each pathway step must carry
const OPT_TOL = 1e-6;

function boundType(glpk, lb, ub) {
  if (lb === ub) return glpk.GLP_FX;
  if (lb <= -BIG && ub >= BIG) return glpk.GLP_FR;
  if (lb <= -BIG) return glpk.GLP_UP;
  if (ub >= BIG) return glpk.GLP_LO;
  return glpk.GLP_DB;
}

// Build the base LP for a GEM + medium. opts:
//   acc            edits key for data.js getEdit
//   extraCols      [{name, stoich:{met:coef}, lb, ub}]  e.g. a product demand
//   stepCons       [{name, vars:[{name,coef}], lb}]     >= lb rows
//   objective      {direction:'max'|'min', vars:[{name,coef}]}
export function buildLP(glpk, gem, mediumBounds, opts = {}) {
  const rows = {};
  for (const m of gem.metabolites)
    rows[m.id] = { name: m.id, vars: [], bnds: { type: glpk.GLP_FX, ub: 0, lb: 0 } };

  const bounds = [];
  for (const r of gem.reactions) {
    const edit = opts.acc ? getEdit(opts.acc, r.id) : null;
    let lb = edit ? edit.lb : r.lb;
    let ub = edit ? edit.ub : r.ub;
    if (r.ex) lb = Object.prototype.hasOwnProperty.call(mediumBounds, r.id) ? mediumBounds[r.id] : 0;
    bounds.push({ name: r.id, type: boundType(glpk, lb, ub), lb, ub });
    for (const [met, coef] of Object.entries(r.stoich))
      if (rows[met]) rows[met].vars.push({ name: r.id, coef });
  }
  for (const c of opts.extraCols || []) {
    bounds.push({ name: c.name, type: boundType(glpk, c.lb, c.ub), lb: c.lb, ub: c.ub });
    for (const [met, coef] of Object.entries(c.stoich))
      if (rows[met]) rows[met].vars.push({ name: c.name, coef });
  }
  const subjectTo = Object.values(rows);
  for (const sc of opts.stepCons || [])
    subjectTo.push({ name: sc.name, vars: sc.vars, bnds: { type: glpk.GLP_LO, lb: sc.lb, ub: 0 } });

  const o = opts.objective || { direction: 'max', vars: [] };
  return {
    name: 'atlas_fba',
    objective: {
      direction: o.direction === 'min' ? glpk.GLP_MIN : glpk.GLP_MAX,
      name: 'obj',
      vars: o.vars,
    },
    subjectTo,
    bounds,
  };
}

async function solve(glpk, lp) {
  const res = await glpk.solve(lp, { msglev: glpk.GLP_MSG_OFF, presol: true });
  const r = res.result;
  return { status: r.status, optimal: r.status === glpk.GLP_OPT, z: r.z, vars: r.vars || {} };
}

// GLPK solution status, in words a reader can act on.
export function statusName(code) {
  return {
    1: 'no feasible solution found',
    2: 'feasible, not proven optimal',
    3: 'infeasible',
    4: 'no feasible solution',
    5: 'optimal',
    6: 'unbounded',
  }[code] || `solver status ${code}`;
}

// Maximum growth (biomass objective) on the medium. Returns {optimal, mu, status}.
export async function maxGrowth(gem, acc, mediumBounds) {
  const glpk = await getGLPK();
  const biomass = gem.stats.biomass_id;
  const lp = buildLP(glpk, gem, mediumBounds, {
    acc, objective: { direction: 'max', vars: [{ name: biomass, coef: 1 }] },
  });
  const s = await solve(glpk, lp);
  return { optimal: s.optimal, mu: s.optimal ? s.z : null, status: s.status };
}

// Product target: the exchange of the _e form when the GEM has one, otherwise a
// demand column on the metabolite itself. Returns {kind, id, extraCols}.
export function productTarget(gem, mid) {
  const base = mid.replace(/_[a-z]+$/, '');
  const exId = `EX_${base}_e`;
  const ex = gem.reactions.find(r => r.ex && r.id === exId);
  if (ex && mid.endsWith('_e')) return { kind: 'exchange', id: exId, extraCols: [] };
  const gemHasMet = gem.metabolites.some(m => m.id === mid);
  if (!gemHasMet) return null;
  return {
    kind: 'demand', id: `DM_${mid}`,
    extraCols: [{ name: `DM_${mid}`, stoich: { [mid]: -1 }, lb: 0, ub: 1000 }],
  };
}

// Directed step constraints for a pathway in this GEM: each step's carried
// alternatives must sum (in the pathway direction) to at least STEP_MIN_FLUX.
// Returns {stepCons, missingStep} where missingStep is the first step index the
// GEM does not carry, or null.
export function stepConstraints(gem, pathway) {
  const gemRxns = new Set(gem.reactions.map(r => r.id));
  const stepCons = [];
  for (let i = 0; i < pathway.steps.length; i++) {
    const vars = [];
    const seen = new Set();
    for (const alt of pathway.steps[i].rxns) {
      if (!gemRxns.has(alt.id) || seen.has(alt.id)) continue;
      seen.add(alt.id);
      vars.push({ name: alt.id, coef: alt.dir === 'rev' ? -1 : 1 });
    }
    if (!vars.length) return { stepCons: null, missingStep: i };
    stepCons.push({ name: `step_${i}`, vars, lb: STEP_MIN_FLUX });
  }
  return { stepCons, missingStep: null };
}

// Feasibility: maximise the product target subject to every pathway step
// carrying flux in the pathway direction. Feasible = optimal LP with product
// flux above tolerance.
export async function pathwayFeasibility(gem, acc, mediumBounds, pathway, target) {
  const glpk = await getGLPK();
  const { stepCons, missingStep } = stepConstraints(gem, pathway);
  if (!stepCons) return { testable: false, missingStep };
  const lp = buildLP(glpk, gem, mediumBounds, {
    acc, extraCols: target.extraCols, stepCons,
    objective: { direction: 'max', vars: [{ name: target.id, coef: 1 }] },
  });
  const s = await solve(glpk, lp);
  const feasible = s.optimal && s.z > OPT_TOL;
  return { testable: true, feasible, productFlux: s.optimal ? s.z : null, status: s.status };
}

// Per-step diagnosis of an infeasible pathway: for each step, maximise the
// product with ONLY that step's constraint active. A step that fails alone
// cannot carry flux in the pathway direction under this GEM's bounds and
// medium; if every step passes alone, the combination is what fails.
export async function diagnoseSteps(gem, acc, mediumBounds, pathway, target) {
  const glpk = await getGLPK();
  const { stepCons } = stepConstraints(gem, pathway);
  if (!stepCons) return null;
  const out = [];
  for (let i = 0; i < stepCons.length; i++) {
    const lp = buildLP(glpk, gem, mediumBounds, {
      acc, extraCols: target.extraCols, stepCons: [stepCons[i]],
      objective: { direction: 'max', vars: [{ name: target.id, coef: 1 }] },
    });
    const s = await solve(glpk, lp);
    out.push({ step: i, ok: s.optimal && s.z > OPT_TOL, flux: s.optimal ? s.z : null });
  }
  return out;
}

// pFBA at the product optimum: fix product flux >= 0.999 * optimum (with the
// step constraints), minimise sum |v|. Returns the flux distribution.
export async function pathwayPFBA(gem, acc, mediumBounds, pathway, target, productOpt) {
  const glpk = await getGLPK();
  const { stepCons } = stepConstraints(gem, pathway);
  if (!stepCons) return { optimal: false };
  const lp = buildLP(glpk, gem, mediumBounds, {
    acc, extraCols: target.extraCols, stepCons,
    objective: { direction: 'min', vars: [] },
  });
  for (const b of lp.bounds) if (b.name === target.id) {
    b.lb = productOpt * 0.999;
    b.ub = Math.max(b.ub, productOpt);
    b.type = boundType(glpk, b.lb, b.ub);
  }
  const absVars = [];
  const extraCons = [];
  for (const r of gem.reactions) {
    const a = 'abs_' + r.id;
    absVars.push({ name: a, coef: 1 });
    lp.bounds.push({ name: a, type: glpk.GLP_LO, lb: 0, ub: BIG });
    extraCons.push({ name: 'ap_' + r.id, vars: [{ name: a, coef: 1 }, { name: r.id, coef: -1 }], bnds: { type: glpk.GLP_LO, lb: 0, ub: 0 } });
    extraCons.push({ name: 'an_' + r.id, vars: [{ name: a, coef: 1 }, { name: r.id, coef: 1 }], bnds: { type: glpk.GLP_LO, lb: 0, ub: 0 } });
  }
  lp.subjectTo = lp.subjectTo.concat(extraCons);
  lp.objective = { direction: glpk.GLP_MIN, name: 'total_flux', vars: absVars };
  const s = await solve(glpk, lp);
  if (!s.optimal) return { optimal: false, status: s.status };
  const fluxes = {};
  for (const r of gem.reactions) fluxes[r.id] = s.vars[r.id] ?? null;
  fluxes[target.id] = s.vars[target.id] ?? null;
  return { optimal: true, fluxes, totalFlux: s.z };
}

// FVA over the given reaction ids with product flux held at >= 0.99 * optimum
// and the step constraints active. onProgress(done, total).
export async function pathwayFVA(gem, acc, mediumBounds, pathway, target, productOpt, rxnIds, onProgress) {
  const glpk = await getGLPK();
  const { stepCons } = stepConstraints(gem, pathway);
  if (!stepCons) return { optimal: false, ranges: {} };
  const lp = buildLP(glpk, gem, mediumBounds, {
    acc, extraCols: target.extraCols, stepCons,
    objective: { direction: 'min', vars: [] },
  });
  for (const b of lp.bounds) if (b.name === target.id) {
    b.lb = productOpt * 0.99;
    b.ub = Math.max(b.ub, productOpt);
    b.type = boundType(glpk, b.lb, b.ub);
  }
  const ranges = {};
  let done = 0;
  for (const rid of rxnIds) {
    lp.objective = { direction: glpk.GLP_MIN, name: 'fva', vars: [{ name: rid, coef: 1 }] };
    const mn = await solve(glpk, lp);
    lp.objective = { direction: glpk.GLP_MAX, name: 'fva', vars: [{ name: rid, coef: 1 }] };
    const mx = await solve(glpk, lp);
    ranges[rid] = {
      min: mn.optimal ? mn.z : null,
      max: mx.optimal ? mx.z : null,
    };
    if (onProgress) onProgress(++done, rxnIds.length);
  }
  return { optimal: true, ranges };
}
