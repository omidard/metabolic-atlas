# Metabolic Atlas — build spec and data contract

A client-side (GitHub Pages) platform to **search bioconversions of a substrate to a product across 35
genome-scale metabolic models (GEMs) of 4 species**, on a 3D metabolic map, with two search modes
(enzyme presence/absence, and flux feasibility), a browsable/editable GEM viewer, and constraint-based
analyses. No server: all compute runs in the browser (JS + GLPK-WASM).

## The one number the landing page must answer
"Which species/strains can convert **<substrate>** to **<product>**, by what pathway, and is it
flux-feasible?" — that answer (a ranked pathway list + the 3D map highlight) is in the first viewport.

## Species / GEMs (35 total)
Parageobacillus thermoglucosidasius (8), Pseudomonas putida (7), Cupriavidus necator (14),
Eubacterium limosum (6). Accents: Parageo `#C0793A`, P.putida `#2E6E8E`, C.necator `#3E8E6E`,
E.limosum `#7A5EA6`.

## Data contract (`docs/data/`)
- **index.json** — `{species:[{key,name,accent,n}], gems:[{acc,species,species_key,accent,genes,reactions,metabolites,n_seqs}]}`
- **media.json** — `{ "<label>": {species, axis, carbon_exchange, carbon_cap, supplements:{ex:lb}, components:{ex:lb}} }` (9 predefined; editable/clonable in the UI)
- **gems/<acc>.json** — one GEM:
  - `stats{genes,reactions,metabolites,exchanges,transporters,mass_balanced,balanceable,ngam,gam,ngam_source,biomass_id}`
  - `growth{<axis>:[muMin,muMax]}` (validation anchors), `biomass{met:coef}`
  - `reactions[{id,name,subsystem,group,lb,ub,gpr,genes[],stoich{met:coef},xr{ec,kegg,rhea,metanetx,bigg},ex,transport}]`
  - `metabolites[{id,name,formula,charge,comp,currency,xr{kegg,chebi,metanetx,seed,inchikey,bigg,hmdb,biocyc}}]`
  - `genes[id]`
- **seqs/<acc>.json** — `{gene_id: nucleotide_sequence}` (lazy-load; ~1-2 MB each)
- **global_graph.json** — the union map across all 35 GEMs:
  - `groups[]` (14 pathway groups), `shells{e:42,p:26,c:13}` (compartment radii)
  - `metabolites{mid:{n(name),c(compartment),g(group),cur(0/1 currency),p:[x,y,z]}}` (3706)
  - `reactions[{id,n,g(group),comp[],sp[species names],s[substrate mids],p[product mids]}]` (5055)

## Tech stack (reuse the institute precedents)
- **3D map**: three.js (WebGL). Render `global_graph.json`: metabolites as nodes on 3 compartment shells
  (outer=exchange/extracellular r=42, middle=periplasm r=26, inner=cytosol r=13); reactions as edges
  substrate→product. Colour by pathway group (one hue per group, a single legend). Currency metabolites
  (`cur:1`) are pre-spread across the shell (positions given) and rendered small/dim so the map is not a
  hairball; offer a "hide currency" toggle. Familiar groups (TCA, PPP, ETC, Glycolysis, Fatty acid,
  Amino acid…) occupy contiguous angular sectors — label each sector. Seed positions are in `p`; a light
  force pass constrained to each metabolite's shell radius may refine, but keep group sectors intact.
  Orbit controls, hover tooltip (name, id, group, xrefs), click to focus.
- **FBA (Mode 2)**: `docs/vendor/glpk.esm.js` (GLPK-WASM, already vendored) — build S·v=0, bounds from
  the GEM + selected medium, objective = biomass or a demand on the product. Same pattern as the
  EcopanGEM Flux Analysis Studio.
- **Pathway enumeration (Mode 1)**: client-side BFS/DFS over the reaction graph (substrate→product),
  same idea as PanRoute (`/data/bioconversion/panroute`). Rank shortest first.

## Search — the core loop
Input: **substrate** + **product** (metabolite pickers with name/BiGG/KEGG search), a species/strain
filter, and a mode.
- **Mode 1 (enzyme presence/absence):** enumerate all pathways substrate→product through reactions,
  ranked by length (shortest = best). For each pathway list the reactions, their enzymes/genes, and
  **which species-strain carries the full pathway** (a presence matrix like PanRoute). Highlight the
  pathway on the 3D map (animate the walk).
- **Mode 2 (flux feasibility):** same, but for a chosen GEM + medium, verify each pathway is flux-
  feasible: set the medium (the substrate REPLACES the carbon source — swap `carbon_exchange` to the
  substrate's exchange), run pFBA/FVA/flux-sampling, mark pathways feasible/infeasible, show the flux
  through each reaction, animate flux on the map.
- **Result panel** (dedicated, right side): ranked pathway cards; each expandable to reactions+genes,
  the species/strain presence matrix, feasibility + flux, and an **Export** (CSV/JSON/SBML).

## Media control
Predefined media selectable; clone + edit (add/remove components, change uptake bounds). When a
substrate is chosen for Mode 2, auto-swap it in as the carbon source (set its exchange lb to the
carbon_cap, close the previous carbon source), and let the user override.

## GEM browser + editor
Load any GEM → show: **stats card** (genes, reactions, metabolites, exchanges, transporters, mass-
balanced/balanceable, GAM, NGAM + source, biomass id), **validation** (growth anchors per axis), the
**biomass reaction formulation**, a searchable **reactions table** (id, name, subsystem, equation,
bounds, GPR, EC, xref links) with **editable lb/ub and directionality** (edits propagate to all
analyses), **genes** (with sequence viewer, from seqs/), and **cross-reference links** (KEGG/ChEBI/
MetaNetX/BiGG/RHEA/SEED). Everything exportable.

## Analyses (constraint-based) — advanced, phase 2
On a GEM + medium: pFBA, FVA, flux sampling; single/double **reaction knockout**; **knockout-to-
feasibility** (find the KO set that makes an infeasible target pathway feasible, or that growth-couples
the product — report shadow prices, biomass change, coupling yield); **MOMA**; over-expression /
down-regulation for growth coupling (advance the FluxStudio idea). All results exportable; all editable
bounds feed back in.

## Visualisation (lots, honest)
3D map (the centrepiece) + per-GEM stat dashboards (bars/donuts with denominators), pathway cards,
presence matrix heatmap, flux-on-map animation, FVA range plots, phenotype phase planes for coupling.
Every count carries its denominator; an absent value renders as "not computed", never a confident 0.

## Design (house rules — non-negotiable)
Light, publication-grade, one accent (per active species), no glow/neon (PanRoute's dark-neon surface
measurably hurt the scientific claim). Operable at 390px; every primary control reachable at every
width; keyboard-operable; WCAG-AA; nothing under 12px. Honest empty/loading/error states. Route through
`.claude/skills/web-design/SKILL.md`; pass that routing into every UI subagent. Verify with
`web_design_gate.py` and `/gate4` before deploy. Name: **Metabolic Atlas**.

## Build phases
1. Data foundation (done) + landing + 3D global map + GEM browser (stats/validation/reactions/genes/
   seqs/xrefs, editable bounds) + Mode-1 search (enumeration + presence matrix + map highlight).
2. Mode-2 flux (GLPK-WASM): medium editor, pFBA/FVA/sampling, feasibility, flux-on-map.
3. Advanced analyses: knockout-to-feasibility, MOMA, growth-coupling, over/under-expression; full export.
