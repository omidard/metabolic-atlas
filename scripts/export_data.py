#!/usr/bin/env python3
"""
Metabolic Atlas - data export. Turns the 35 GEMs_p9 into the web platform's data contract.
Outputs to docs/data/: index.json, media.json, gems/<acc>.json, seqs/<acc>.json, global_graph.json.
"""
import cobra, json, glob, os, re, warnings, math
warnings.filterwarnings("ignore")
V2="/data/metabolic_atlas_v2"
RX="/data/Brilliant_genomics_department/brilliant_software_dev/Reactome/data/generated"
OUT="/data/metabolic_atlas_platform/docs/data"
for d in ("gems","seqs"): os.makedirs(f"{OUT}/{d}",exist_ok=True)

SPECIES={"Parageobacillus":("Parageobacillus thermoglucosidasius","#C0793A"),
         "Pseudomonas_putida":("Pseudomonas putida","#2E6E8E"),
         "Cupriavidus_necator":("Cupriavidus necator","#3E8E6E"),
         "Eubacterium_limosum":("Eubacterium limosum","#7A5EA6")}
NGAM={"Parageobacillus":(3.141,152.3,"MOL2021 fitted"),"Pseudomonas_putida":(3.4,46.9,"van Duuren 2013 / iJN1463"),
      "Cupriavidus_necator":(3.0,682.4,"iCN1361"),"Eubacterium_limosum":(2.0,195.1,"iHN637-anchored")}
# validation growth anchors (from the recuration)
GROWTH={"Parageobacillus":{"glucose":[0.629,0.629],"acetate":[0.39,0.391]},
        "Pseudomonas_putida":{"glucose":[0.353,1.085],"acetate":[0.344,0.958]},
        "Cupriavidus_necator":{"fructose":[0.248,0.274],"autotrophic":[0.895,0.901]},
        "Eubacterium_limosum":{"glucose":[0.233,0.233],"methanol":[0.100,0.100]}}
CURRENCY={"h","h2o","atp","adp","amp","pi","ppi","nad","nadh","nadp","nadph","co2","o2","coa",
          "nh4","h2","fad","fadh2","q8","q8h2","mqn8","mql8","fdxo_42","fdxrd_42","gtp","gdp","gmp",
          "utp","udp","ump","ctp","cdp","cmp","itp","idp","so4","so3","h2s","thf","mlthf","10fthf",
          "amet","ahcys","acp","actp","glu__L","akg","na1","k","cl","mg2","fe2","fe3","mn2","zn2",
          "cu2","ca2","cobalt2","ni2","h2o2","pyr","accoa"}
SUBSYS_GROUP={  # BiGG subsystem -> display pathway group (familiar patterns)
  "Citric Acid Cycle":"TCA","Glycolysis/Gluconeogenesis":"Glycolysis","Pentose Phosphate Pathway":"PPP",
  "Oxidative Phosphorylation":"ETC","Fatty Acid Biosynthesis":"Fatty acid","Fatty Acid Synthesis":"Fatty acid",
  "Purine and Pyrimidine Biosynthesis":"Nucleotide","Nucleotide Salvage Pathway":"Nucleotide",
  "Cofactor and Prosthetic Group Biosynthesis":"Cofactor","Cell Envelope Biosynthesis":"Cell envelope",
  "Membrane Lipid Metabolism":"Lipid","Glycerophospholipid Metabolism":"Lipid","Transport, Inner Membrane":"Transport",
  "Transport, Outer Membrane Porin":"Transport","Transport, Outer Membrane":"Transport",
  "Alternate Carbon Metabolism":"Carbon","Inorganic Ion Transport and Metabolism":"Ion transport"}
AA={"Alanine","Arginine","Asparagine","Aspartate","Cysteine","Glutamate","Glutamine","Glycine","Histidine",
    "Isoleucine","Leucine","Lysine","Methionine","Phenylalanine","Proline","Serine","Threonine","Tryptophan","Tyrosine","Valine"}

# reaction -> subsystem (consensus across BiGG models)
def load_subsys():
    from collections import Counter
    c=Counter(); sub=defaultdict(Counter)
    with open(f"{RX}/model_reactions_prok.tsv") as f:
        next(f)
        for line in f:
            p=line.rstrip("\n").split("\t")
            if len(p)>=5 and p[4]: sub[p[0]][p[4]]+=1
    return {r:cc.most_common(1)[0][0] for r,cc in sub.items()}
from collections import defaultdict, Counter
SUBSYS=load_subsys()
def subgroup(sname):
    if not sname: return None
    if sname in SUBSYS_GROUP: return SUBSYS_GROUP[sname]
    for aa in AA:
        if aa.lower() in sname.lower(): return "Amino acid"
    if "amino acid" in sname.lower(): return "Amino acid"
    return None

def xrefs_r(ann):
    keep={"ec-code":"ec","kegg.reaction":"kegg","rhea":"rhea","metanetx.reaction":"metanetx","bigg.reaction":"bigg"}
    o={}
    for k,v in (ann or {}).items():
        if k in keep: o[keep[k]]=(v if isinstance(v,list) else [v])
    return o
def xrefs_m(ann):
    keep={"kegg.compound":"kegg","chebi":"chebi","metanetx.chemical":"metanetx","seed.compound":"seed",
          "inchi_key":"inchikey","bigg.metabolite":"bigg","hmdb":"hmdb","biocyc":"biocyc"}
    o={}
    for k,v in (ann or {}).items():
        if k in keep: o[keep[k]]=(v if isinstance(v,list) else [v])
    return o

def export_gem(acc, sp):
    m=cobra.io.read_sbml_model(f"{V2}/{sp}/GEMs_p9/{acc}.xml")
    bm=[r for r in m.reactions if 'BIOMASS' in r.id.upper()][0]
    def base(rid): return rid.replace("R_","")
    rxns=[]
    for r in m.reactions:
        sub=SUBSYS.get(base(r.id),"")
        rxns.append(dict(id=r.id, name=r.name or r.id, subsystem=sub, group=subgroup(sub),
                         lb=r.lower_bound, ub=r.upper_bound, gpr=r.gene_reaction_rule,
                         genes=[g.id for g in r.genes],
                         stoich={met.id:coef for met,coef in r.metabolites.items()},
                         xr=xrefs_r(r.annotation),
                         ex=r.id.startswith(("EX_","R_EX_")),
                         transport=(len({mm.id.rsplit("_",1)[-1] for mm in r.metabolites})>1 and not r.id.startswith(("EX_","R_EX_")))))
    mets=[]
    for mt in m.metabolites:
        mets.append(dict(id=mt.id, name=mt.name or mt.id, formula=mt.formula or "", charge=mt.charge or 0,
                         comp=mt.compartment or (mt.id.rsplit("_",1)[-1] if "_" in mt.id else "c"),
                         currency=(re.sub(r'_[a-z0-9]+$','',mt.id) in CURRENCY),
                         xr=xrefs_m(mt.annotation)))
    n_ex=sum(1 for r in m.reactions if r.id.startswith(("EX_","R_EX_")))
    balanceable=[r for r in m.reactions if not r.id.startswith(("EX_","R_EX_","DM_","SK_")) and 'BIOMASS' not in r.id.upper()
                 and all(mm.formula for mm in r.metabolites)]
    bal=sum(1 for r in balanceable if not r.check_mass_balance())
    ngam,gam,ngsrc=NGAM[sp]
    data=dict(acc=acc, species=SPECIES[sp][0], accent=SPECIES[sp][1],
              stats=dict(genes=len(m.genes), reactions=len(m.reactions), metabolites=len(m.metabolites),
                         exchanges=n_ex, transporters=sum(1 for r in rxns if r["transport"]),
                         mass_balanced=bal, balanceable=len(balanceable),
                         ngam=ngam, gam=gam, ngam_source=ngsrc, biomass_id=bm.id),
              biomass={mt.id:round(coef,6) for mt,coef in bm.metabolites.items()},
              growth=GROWTH[sp], reactions=rxns, metabolites=mets, genes=[g.id for g in m.genes])
    json.dump(data, open(f"{OUT}/gems/{acc}.json","w"), separators=(",",":"))
    # sequences from .ffn (header id -> seq)
    seqs={}
    ffn=f"{V2}/{sp}/input_annotations/{acc}/{acc}.ffn"
    if os.path.exists(ffn):
        cur=None; buf=[]
        model_genes=set(g.id for g in m.genes)
        for line in open(ffn):
            if line.startswith(">"):
                if cur and cur in model_genes: seqs[cur]="".join(buf)
                cur=line[1:].split()[0].strip(); buf=[]
            else: buf.append(line.strip())
        if cur and cur in model_genes: seqs[cur]="".join(buf)
    json.dump(seqs, open(f"{OUT}/seqs/{acc}.json","w"), separators=(",",":"))
    return data["stats"], len(seqs)

if __name__=="__main__":
    index={"species":[],"gems":[]}
    for sp,(disp,accent) in SPECIES.items():
        accs=sorted(os.path.basename(p)[:-4] for p in glob.glob(f"{V2}/{sp}/GEMs_p9/*.xml"))
        index["species"].append(dict(key=sp,name=disp,accent=accent,n=len(accs)))
        for acc in accs:
            st,nseq=export_gem(acc,sp)
            index["gems"].append(dict(acc=acc,species=disp,species_key=sp,accent=accent,
                                      genes=st["genes"],reactions=st["reactions"],metabolites=st["metabolites"],n_seqs=nseq))
            print(f"  {acc} ({disp}): {st['reactions']} rxn, {st['genes']} genes, {nseq} seqs")
    json.dump(index, open(f"{OUT}/index.json","w"), indent=1)
    print(f"\nexported {len(index['gems'])} GEMs")
