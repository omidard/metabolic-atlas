#!/usr/bin/env python3
"""Build docs/assets/graph_meta.json from docs/data/ (read-only inputs).

Derived sidecar for the client app:
  - accs[]: the 35 GEM accessions in fixed order (bit order of the presence mask)
  - rxns{id}: m = hex presence bitmask over accs (bit i set = accs[i] carries the
    reaction id), d = direction union across carrying GEMs (1 fwd only, 2 rev only,
    3 both), e = union of EC numbers from GEM cross-references
Re-run whenever docs/data/ changes. Inputs are never modified.
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "docs", "data")
OUT = os.path.join(ROOT, "docs", "assets", "graph_meta.json")


def main():
    index = json.load(open(os.path.join(DATA, "index.json")))
    graph = json.load(open(os.path.join(DATA, "global_graph.json")))
    graph_ids = [r["id"] for r in graph["reactions"]]
    accs = [
        {"acc": g["acc"], "sk": g["species_key"], "sp": g["species"]}
        for g in index["gems"]
    ]

    rxns = {rid: {"mask": 0, "fwd": False, "rev": False, "ec": set()} for rid in graph_ids}
    unmatched_gem_rxns = 0
    for i, a in enumerate(accs):
        gem = json.load(open(os.path.join(DATA, "gems", a["acc"] + ".json")))
        for r in gem["reactions"]:
            rec = rxns.get(r["id"])
            if rec is None:
                unmatched_gem_rxns += 1
                continue
            rec["mask"] |= 1 << i
            if r["ub"] > 0:
                rec["fwd"] = True
            if r["lb"] < 0:
                rec["rev"] = True
            for ec in r.get("xr", {}).get("ec", []) or []:
                rec["ec"].add(ec)

    out_rxns = {}
    orphan_graph_rxns = 0
    for rid, rec in rxns.items():
        if rec["mask"] == 0:
            orphan_graph_rxns += 1
        d = (1 if rec["fwd"] else 0) | (2 if rec["rev"] else 0)
        if d == 0:
            d = 1  # no carrying GEM found; traverse as written
        entry = {"m": format(rec["mask"], "x"), "d": d}
        if rec["ec"]:
            entry["e"] = sorted(rec["ec"])
        out_rxns[rid] = entry

    payload = {
        "generated_from": "docs/data/index.json + docs/data/gems/*.json",
        "n_accs": len(accs),
        "n_rxns": len(out_rxns),
        "accs": accs,
        "rxns": out_rxns,
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    print(f"wrote {OUT} ({os.path.getsize(OUT)/1024:.0f} KiB)")
    print(f"graph reactions: {len(graph_ids)}; carried by no GEM: {orphan_graph_rxns}")
    print(f"GEM reactions not in graph (skipped): {unmatched_gem_rxns}")
    rev = sum(1 for r in out_rxns.values() if r["d"] & 2)
    print(f"reactions with a reverse direction in at least one GEM: {rev}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
