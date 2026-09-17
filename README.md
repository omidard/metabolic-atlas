# Metabolic Atlas

A client-side platform to **search bioconversions of a substrate to a product across 35 genome-scale
metabolic models (GEMs) of 4 species**, on a 3D metabolic map, with enzyme-presence and flux-feasibility
search modes, a browsable/editable GEM viewer, constraint-based analyses, sequences and cross-references.

- **Species (35 GEMs):** Parageobacillus thermoglucosidasius (8), Pseudomonas putida (7),
  Cupriavidus necator (14), Eubacterium limosum (6).
- **Runs entirely in the browser** (GitHub Pages; GLPK-WASM for flux analysis; three.js for the 3D map).
- Data foundation and the build spec: see `BUILD_SPEC.md`. Data lives in `docs/data/`.

Not a KEGG-map tool: the map is a global 3D network built only from these species' reactions, layered by
compartment (exchange/periplasm/cytosol) with familiar pathway sectors (TCA, PPP, ETC, glycolysis,
fatty-acid and amino-acid biosynthesis).

## Data (docs/data/)
`index.json` (GEM index) · `media.json` (9 editable predefined media) · `gems/<acc>.json` (per GEM:
stats, validation, reactions, metabolites, GPRs, bounds, cross-refs) · `seqs/<acc>.json` (gene
sequences) · `global_graph.json` (the 3D map: 3,706 metabolites, 5,055 reactions, 14 pathway groups).

Regenerate: `python scripts/export_data.py && python scripts/export_media_graph.py`.
