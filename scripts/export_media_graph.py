#!/usr/bin/env python3
"""Export media.json (editable predefined media) + global_graph.json (3D layered metabolic map)."""
import cobra, json, glob, os, re, warnings, math
warnings.filterwarnings("ignore")
from collections import defaultdict, Counter
V2="/data/metabolic_atlas_v2"; V1="/data/metabolic_atlas"; OUT="/data/metabolic_atlas_platform/docs/data"
import sys; sys.path.insert(0,f"{V2}/scripts")
from atlas_species_config import SPECIES_CONFIG

# ---- media.json: the validation media, as editable {exchange: uptake_bound} ----
media={}
for sp,cfg in SPECIES_CONFIG.items():
    for axis,(mf,cx,cap) in cfg["media"].items():
        try: comp=json.load(open(f"{cfg['media_dir']}/{mf}"))
        except Exception: continue
        key=f"{sp.split('_')[0]} - {axis}"
        media[key]=dict(species=sp, axis=axis, carbon_exchange=cx, carbon_cap=cap,
                        supplements=cfg["supplements"], components=comp)
json.dump(media, open(f"{OUT}/media.json","w"), indent=1)
print(f"media.json: {len(media)} predefined media")

# ---- global_graph.json: union of reactions across all 35 GEMs, 3D layered ----
CURRENCY={"h","h2o","atp","adp","amp","pi","ppi","nad","nadh","nadp","nadph","co2","o2","coa","nh4","h2",
          "fad","fadh2","q8","q8h2","mqn8","mql8","gtp","gdp","utp","udp","ctp","cdp","so4","thf","amet","ahcys"}
GROUP_ORDER=["Glycolysis","TCA","PPP","ETC","Fatty acid","Amino acid","Nucleotide","Cofactor","Lipid",
             "Cell envelope","Carbon","Ion transport","Transport","Other"]
# reuse subsystem map from the exported gems (group already computed per reaction)
rxn_group={}; rxn_name={}; rxn_comp=defaultdict(set); rxn_species=defaultdict(set); rxn_stoich={}
met_name={}; met_comp={}; met_currency={}
for gp in glob.glob(f"{OUT}/gems/*.json"):
    d=json.load(open(gp)); spk=d["species"]
    for r in d["reactions"]:
        b=r["id"].replace("R_","")
        rxn_group.setdefault(b, r.get("group") or "Other"); rxn_name.setdefault(b,r["name"])
        rxn_species[b].add(spk)
        if b not in rxn_stoich: rxn_stoich[b]=r["stoich"]
        for mid in r["stoich"]:
            comp=mid.rsplit("_",1)[-1]
            if r.get("ex"): rxn_comp[b].add("e")
            else: rxn_comp[b].add(comp)
    for mt in d["metabolites"]:
        met_name.setdefault(mt["id"],mt["name"]); met_comp.setdefault(mt["id"],mt["comp"])
        met_currency[mt["id"]]=mt["currency"]

# dominant group per metabolite (by the reactions it participates in)
met_group=defaultdict(Counter)
for b,st in rxn_stoich.items():
    g=rxn_group.get(b,"Other")
    for mid in st: met_group[mid][g]+=1
def dgroup(mid):
    c=met_group.get(mid)
    return c.most_common(1)[0][0] if c else "Other"

# 3D layout: compartment shell (radius) x subsystem angular sector
SHELL={"e":42.0,"p":26.0,"c":13.0}
groups=[g for g in GROUP_ORDER if any(dgroup(m)==g for m in met_name)]
sector={g:i for i,g in enumerate(groups)}
NS=max(len(groups),1)
def pos(mid,i):
    comp=met_comp.get(mid,"c"); comp = comp if comp in SHELL else "c"
    r=SHELL[comp]
    g=dgroup(mid); s=sector.get(g,0)
    # currency metabolites: spread evenly around the whole shell (not clustered)
    if met_currency.get(mid):
        phi=(hash(mid)%997)/997.0*2*math.pi; th=(hash(mid+"z")%997)/997.0*math.pi
    else:
        base=(s+0.5)/NS*2*math.pi
        phi=base + ((hash(mid)%200)/200.0-0.5)*(2*math.pi/NS)*0.85
        th=math.pi*(0.28 + 0.44*((hash(mid+"t")%200)/200.0))
    x=r*math.sin(th)*math.cos(phi); y=r*math.cos(th); z=r*math.sin(th)*math.sin(phi)
    return [round(x,2),round(y,2),round(z,2)]

# emit compact: metabolites (non-currency get positions; currency flagged), reactions (edges)
mnodes={}; i=0
for mid in sorted(met_name):
    mnodes[mid]=dict(n=met_name[mid], c=met_comp.get(mid,"c"), g=dgroup(mid),
                     cur=1 if met_currency.get(mid) else 0, p=pos(mid,i)); i+=1
redges=[]
for b in sorted(rxn_stoich):
    st=rxn_stoich[b]
    subs=[m for m,c in st.items() if c<0]; prods=[m for m,c in st.items() if c>0]
    redges.append(dict(id=b, n=rxn_name.get(b,b), g=rxn_group.get(b,"Other"),
                       comp=sorted(rxn_comp.get(b,["c"])), sp=sorted(rxn_species.get(b,[])),
                       s=subs, p=prods))
graph=dict(groups=groups, shells=SHELL, metabolites=mnodes, reactions=redges,
           n_metabolites=len(mnodes), n_reactions=len(redges))
json.dump(graph, open(f"{OUT}/global_graph.json","w"), separators=(",",":"))
print(f"global_graph.json: {len(mnodes)} metabolite nodes, {len(redges)} reactions, {len(groups)} pathway groups")
print("groups:", groups)
