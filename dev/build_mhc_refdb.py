"""Reference sequences for 'closest allele / gene' of MHC chains (imm_annotate).

Sources: IPD-IMGT/HLA (all proteins, 2-field names), IPD-MHC (non-human alleles), UniProt
MHC class I / II / B2M entries (incl. InterPro class II).  For every sequence the
peptide-binding region (class I alpha1+alpha2, class II alpha1 or beta1) or the beta2m
domain is cut out with the MHC HMMs; identical regions are merged (names joined).

    python dev/build_mhc_refdb.py SEQDIR immtools/models/mhc_refs.fasta.gz
"""
import collections, gzip, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from immtools import _mhc as H
import seqsources as C

REGION = {"MHC_I_ALPHA": ("G_ALPHA1", "G_ALPHA2"), "MHC_II_ALPHA": ("G_ALPHA", "G_ALPHA"),
          "MHC_II_BETA": ("G_BETA", "G_BETA"), "B2M": ("C_B2M", "C_B2M")}
IPD_SPECIES = {
    "Mafa": "Macaca fascicularis", "Mamu": "Macaca mulatta", "Mane": "Macaca nemestrina",
    "BoLA": "Bos taurus", "DLA": "Canis lupus familiaris", "SLA": "Sus scrofa",
    "Patr": "Pan troglodytes", "Paan": "Papio anubis", "Ovar": "Ovis aries",
    "Caja": "Callithrix jacchus", "Chsa": "Chlorocebus sabaeus", "Aona": "Aotus nancymaae",
    "Ceat": "Cercocebus atys", "Math": "Macaca thibetana", "Chpy": "Chlorocebus pygerythrus",
    "Gogo": "Gorilla gorilla", "Malo": "Macaca leonina", "Sasa": "Salmo salar",
    "Eqca": "Equus caballus", "Rano": "Rattus norvegicus", "Popy": "Pongo pygmaeus",
    "Papa": "Pan paniscus", "Saoe": "Saguinus oedipus", "Onmy": "Oncorhynchus mykiss",
    "Poab": "Pongo abelii", "Chae": "Chlorocebus aethiops", "Bogr": "Bos grunniens",
    "Bubu": "Bubalus bubalis", "Bofr": "Bos frontalis", "Aovo": "Aotus vociferans",
    "Aoni": "Aotus nigriceps", "Maar": "Macaca arctoides", "Cahi": "Capra hircus",
    "Mafu": "Macaca fuscata", "Gaga": "Gallus gallus", "Paur": "Papio ursinus",
    "Atfu": "Ateles fusciceps", "Maas": "Macaca assamensis", "Ctid": "Ctenopharyngodon idella",
    "Eqas": "Equus asinus", "Feca": "Felis catus", "Susc": "Sus scrofa"}
PRIORITY = {"IPD-IMGT/HLA": 0, "IPD-MHC": 1, "UniProt": 2}


def hla_all(seqdir):
    seen = set()
    for name, seq in C.read_fasta(seqdir / "hla_prot.fasta"):
        m = re.match(r"HLA:(\S+)\s+(\S+)", name)
        if not m:
            continue
        allele = "HLA-" + ":".join(m.group(2).split(":")[:2])
        if (allele, seq) in seen:
            continue
        seen.add((allele, seq))
        yield dict(uid=f"HLA:{m.group(1)}", source="IPD-IMGT/HLA", acc=m.group(1),
                   species="Homo sapiens", name=allele, seq=seq)


def main(seqdir, out):
    seqdir = Path(seqdir)
    recs = [r for r in C.uniprot(seqdir)] + list(C.ipd_mhc(seqdir)) + list(hla_all(seqdir))
    for r in recs:
        if r["source"] == "IPD-MHC":
            r["species"] = IPD_SPECIES.get(r["species"], r["species"] + " (IPD-MHC code)")
        if r["source"] == "UniProt":
            r["name"] = f"{r['name'].split('(')[0].strip()} [UniProt {r['acc']}]"
    print(len(recs), "sequences", file=sys.stderr)
    clean, hits = C.batch_domains(recs)
    region = collections.OrderedDict()
    for i, r in enumerate(recs):
        hs = hits.get(i, [])
        best = {}
        for h in hs:
            if h["model"] not in best or h["score"] > best[h["model"]]["score"]:
                best[h["model"]] = h
        for ct, (a, b) in REGION.items():
            if a in best and b in best and best[a]["start"] <= best[b]["end"]:
                s = clean[i][best[a]["start"]:best[b]["end"] + 1]
                if len(s) < 60:
                    continue
                key = (ct, s)
                region.setdefault(key, []).append((PRIORITY[r["source"]], r["name"], r["species"]))
                break
    with gzip.open(out, "wt") as fh:
        for (ct, s), names in region.items():
            names.sort()
            sp = collections.Counter(n[2] for n in names).most_common(1)[0][0]
            shown = list(dict.fromkeys(n[1] for n in names))
            label = " / ".join(shown[:3]) + (f" (+{len(shown) - 3})" if len(shown) > 3 else "")
            fh.write(f">{ct}|{label}|{sp}\n{s}\n")
    print(collections.Counter(k[0] for k in region), file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:3])
