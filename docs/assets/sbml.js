// SBML (level 3 version 1 + fbc version 2) and COBRA-JSON export of one GEM,
// with the session's in-memory bound edits applied. Pure string/object
// builders, no DOM, so the node harness can verify the output against the
// source GEM. Absent values are omitted, never invented: a metabolite without
// a strictly elemental formula gets no fbc:chemicalFormula, a reaction without
// a gene rule gets no geneProductAssociation.

import { getEdit, listEdits } from './data.js';

// SBML SId: [A-Za-z_][A-Za-z0-9_]*
const sid = (s) => {
  let t = String(s).replace(/[^A-Za-z0-9_]/g, '_');
  if (!/^[A-Za-z_]/.test(t)) t = '_' + t;
  return t;
};
const xesc = (s) => String(s ?? '').replace(/[&<>"']/g, c => (
  { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
// fbc:chemicalFormula accepts only element-count runs; averaged or dotted
// formulas (C38.14..., C9H23N3O6.Fe) are omitted from the SBML rather than
// written as invalid values. The COBRA JSON keeps the raw string.
const FORMULA_RE = /^([A-Z][a-z]?[0-9]*)+$/;

const COMP_NAMES = {
  c: 'cytosol', e: 'extracellular space', p: 'periplasm', m: 'mitochondrion',
  n: 'nucleus', x: 'peroxisome', f: 'flagellum', h: 'chloroplast',
};

export function effectiveBounds(gem, acc) {
  const out = new Map();
  for (const r of gem.reactions) {
    const e = acc ? getEdit(acc, r.id) : null;
    out.set(r.id, e ? { lb: e.lb, ub: e.ub, edited: true } : { lb: r.lb, ub: r.ub, edited: false });
  }
  return out;
}

// ---- GPR boolean parser: expr := term ('or' term)*; term := fac ('and' fac)*;
//      fac := '(' expr ')' | gene. Returns {gene} | {op:'and'|'or', kids} | null.
export function parseGPR(str) {
  if (!str || !str.trim()) return null;
  const toks = str.replace(/\(/g, ' ( ').replace(/\)/g, ' ) ').trim().split(/\s+/);
  let i = 0;
  const peek = () => toks[i];
  const eat = () => toks[i++];
  function expr() {
    const kids = [term()];
    while (peek() && peek().toLowerCase() === 'or') { eat(); kids.push(term()); }
    return kids.length === 1 ? kids[0] : { op: 'or', kids };
  }
  function term() {
    const kids = [fac()];
    while (peek() && peek().toLowerCase() === 'and') { eat(); kids.push(fac()); }
    return kids.length === 1 ? kids[0] : { op: 'and', kids };
  }
  function fac() {
    if (peek() === '(') {
      eat();
      const e = expr();
      if (eat() !== ')') throw new Error('unbalanced parenthesis');
      return e;
    }
    const t = eat();
    if (t === undefined || t === ')' || /^(and|or)$/i.test(t)) throw new Error(`unexpected token "${t}"`);
    return { gene: t };
  }
  try {
    const tree = expr();
    if (i !== toks.length) throw new Error('trailing tokens');
    return tree;
  } catch {
    return null;   // unparseable rule: the caller counts and reports it
  }
}

function gprGenes(node, set) {
  if (!node) return set;
  if (node.gene) { set.add(node.gene); return set; }
  for (const k of node.kids) gprGenes(k, set);
  return set;
}

function gprXML(node, ind) {
  if (node.gene) return `${ind}<fbc:geneProductRef fbc:geneProduct="G_${sid(node.gene)}"/>`;
  const tag = node.op === 'and' ? 'fbc:and' : 'fbc:or';
  return `${ind}<${tag}>\n${node.kids.map(k => gprXML(k, ind + '  ')).join('\n')}\n${ind}</${tag}>`;
}

// ---- SBML builder. Returns {xml, nEdits, nGprSkipped, nFormulaSkipped}.
export function buildSBML(gem, acc, meta = {}) {
  const bounds = effectiveBounds(gem, acc);
  const edits = acc ? listEdits(acc) : [];
  const modelId = 'M_' + sid(acc || 'model');

  // compartments actually used
  const comps = new Map();
  for (const m of gem.metabolites) {
    const c = m.comp || 'c';
    if (!comps.has(c)) comps.set(c, COMP_NAMES[c.replace(/^C_/, '')] || c);
  }

  // shared bound parameters, one per distinct value
  const paramId = new Map();
  const pFor = (v) => {
    if (!paramId.has(v)) paramId.set(v, `B${paramId.size}`);
    return paramId.get(v);
  };
  for (const b of bounds.values()) { pFor(b.lb); pFor(b.ub); }

  // gene products: the model's gene list plus any id appearing in a rule
  const geneSet = new Set(gem.genes);
  let nGprSkipped = 0;
  const gprTrees = new Map();
  for (const r of gem.reactions) {
    if (r.gpr && r.gpr.trim()) {
      const t = parseGPR(r.gpr);
      if (t) { gprTrees.set(r.id, t); gprGenes(t, geneSet); }
      else nGprSkipped++;
    }
  }

  let nFormulaSkipped = 0;
  const speciesXML = gem.metabolites.map(m => {
    const fm = (m.formula || '').trim();
    const okF = fm && FORMULA_RE.test(fm);
    if (fm && !okF) nFormulaSkipped++;
    const attrs = [
      `id="M_${sid(m.id)}"`,
      m.name ? `name="${xesc(m.name)}"` : '',
      `compartment="${sid(m.comp || 'c')}"`,
      'hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"',
      Number.isInteger(m.charge) ? `fbc:charge="${m.charge}"` : '',
      okF ? `fbc:chemicalFormula="${fm}"` : '',
    ].filter(Boolean).join(' ');
    return `      <species ${attrs}/>`;
  }).join('\n');

  const rxnXML = gem.reactions.map(r => {
    const b = bounds.get(r.id);
    const reac = [], prod = [];
    for (const [mid, coef] of Object.entries(r.stoich)) {
      const line = `          <speciesReference species="M_${sid(mid)}" stoichiometry="${Math.abs(coef)}" constant="true"/>`;
      (coef < 0 ? reac : prod).push(line);
    }
    const tree = gprTrees.get(r.id);
    return [
      `      <reaction id="R_${sid(r.id)}"${r.name ? ` name="${xesc(r.name)}"` : ''} reversible="${b.lb < 0}" fast="false" fbc:lowerFluxBound="${pFor(b.lb)}" fbc:upperFluxBound="${pFor(b.ub)}">`,
      reac.length ? `        <listOfReactants>\n${reac.join('\n')}\n        </listOfReactants>` : '',
      prod.length ? `        <listOfProducts>\n${prod.join('\n')}\n        </listOfProducts>` : '',
      tree ? `        <fbc:geneProductAssociation>\n${gprXML(tree, '          ')}\n        </fbc:geneProductAssociation>` : '',
      '      </reaction>',
    ].filter(Boolean).join('\n');
  }).join('\n');

  const notes = [
    `Exported from Metabolic Atlas${meta.release ? `, data release ${meta.release}` : ''}.`,
    `Species: ${gem.species || 'not recorded'}. Accession: ${acc}.`,
    edits.length
      ? `Session bound edits applied to ${edits.length} of ${gem.reactions.length} reactions: ${edits.map(e => `${e.rid} [${e.lb}, ${e.ub}]`).join('; ')}.`
      : `No session bound edits; bounds are the model's own.`,
    nGprSkipped ? `${nGprSkipped} of ${gem.reactions.length} gene rules could not be parsed and are omitted from the SBML (the source rule strings remain in the JSON export).` : '',
  ].filter(Boolean);

  const xml = `<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2" level="3" version="1" fbc:required="false">
  <model id="${modelId}"${gem.species ? ` name="${xesc(gem.species + ' ' + acc)}"` : ''} fbc:strict="true">
    <notes>
      <body xmlns="http://www.w3.org/1999/xhtml">
${notes.map(n => `        <p>${xesc(n)}</p>`).join('\n')}
      </body>
    </notes>
    <listOfUnitDefinitions>
      <unitDefinition id="mmol_per_gDW_per_hr">
        <listOfUnits>
          <unit kind="mole" exponent="1" scale="-3" multiplier="1"/>
          <unit kind="gram" exponent="-1" scale="0" multiplier="1"/>
          <unit kind="second" exponent="-1" scale="0" multiplier="3600"/>
        </listOfUnits>
      </unitDefinition>
    </listOfUnitDefinitions>
    <listOfCompartments>
${[...comps].map(([c, nm]) => `      <compartment id="${sid(c)}" name="${xesc(nm)}" constant="true"/>`).join('\n')}
    </listOfCompartments>
    <listOfSpecies>
${speciesXML}
    </listOfSpecies>
    <listOfParameters>
${[...paramId].map(([v, id]) => `      <parameter id="${id}" value="${v}" units="mmol_per_gDW_per_hr" constant="true" sboTerm="SBO:0000625"/>`).join('\n')}
    </listOfParameters>
    <listOfReactions>
${rxnXML}
    </listOfReactions>
    <fbc:listOfObjectives fbc:activeObjective="obj">
      <fbc:objective fbc:id="obj" fbc:type="maximize">
        <fbc:listOfFluxObjectives>
          <fbc:fluxObjective fbc:reaction="R_${sid(gem.stats.biomass_id)}" fbc:coefficient="1"/>
        </fbc:listOfFluxObjectives>
      </fbc:objective>
    </fbc:listOfObjectives>
    <fbc:listOfGeneProducts>
${[...geneSet].map(g => `      <fbc:geneProduct fbc:id="G_${sid(g)}" fbc:label="${xesc(g)}"/>`).join('\n')}
    </fbc:listOfGeneProducts>
  </model>
</sbml>
`;
  return { xml, nEdits: edits.length, nGprSkipped, nFormulaSkipped };
}

// ---- COBRA-style JSON (cobrapy schema version 1). Ids stay as in the source
// GEM; xrefs become identifier-style annotations.
const MET_XR_KEYS = {
  kegg: 'kegg.compound', chebi: 'chebi', metanetx: 'metanetx.chemical',
  seed: 'seed.compound', inchikey: 'inchi_key', bigg: 'bigg.metabolite',
  hmdb: 'hmdb', biocyc: 'biocyc',
};
const RXN_XR_KEYS = {
  ec: 'ec-code', kegg: 'kegg.reaction', rhea: 'rhea',
  metanetx: 'metanetx.reaction', bigg: 'bigg.reaction',
};

function xrAnnotation(xr, keyMap) {
  const out = {};
  for (const [k, ids] of Object.entries(xr || {})) {
    if (keyMap[k] && ids && ids.length) out[keyMap[k]] = ids;
  }
  return out;
}

export function buildCobraJSON(gem, acc, meta = {}) {
  const bounds = effectiveBounds(gem, acc);
  const edits = acc ? listEdits(acc) : [];
  const obj = {
    id: sid(acc || 'model'),
    name: `${gem.species || ''} ${acc || ''}`.trim(),
    version: '1',
    compartments: {},
    metabolites: gem.metabolites.map(m => ({
      id: m.id,
      name: m.name || '',
      compartment: m.comp || 'c',
      charge: Number.isInteger(m.charge) ? m.charge : null,
      formula: (m.formula || '').trim() || null,
      notes: {},
      annotation: xrAnnotation(m.xr, MET_XR_KEYS),
    })),
    reactions: gem.reactions.map(r => {
      const b = bounds.get(r.id);
      return {
        id: r.id,
        name: r.name || '',
        metabolites: r.stoich,
        lower_bound: b.lb,
        upper_bound: b.ub,
        gene_reaction_rule: r.gpr || '',
        subsystem: r.subsystem || '',
        objective_coefficient: r.id === gem.stats.biomass_id ? 1 : 0,
        notes: b.edited ? { atlas_session_edit: true } : {},
        annotation: xrAnnotation(r.xr, RXN_XR_KEYS),
      };
    }),
    genes: gem.genes.map(g => ({ id: g, name: g, notes: {}, annotation: {} })),
    notes: {
      exported_from: `Metabolic Atlas${meta.release ? ` data release ${meta.release}` : ''}`,
      species: gem.species || null,
      accession: acc,
      session_bound_edits: edits.length ? edits : [],
    },
  };
  for (const m of gem.metabolites) {
    const c = m.comp || 'c';
    if (!(c in obj.compartments)) obj.compartments[c] = COMP_NAMES[c.replace(/^C_/, '')] || c;
  }
  return { obj, nEdits: edits.length };
}
