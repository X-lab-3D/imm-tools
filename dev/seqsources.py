# SPDX-License-Identifier: Apache-2.0
"""Readers for the MHC sequence sources (UniProt TSV exports, IPD-MHC and IPD-IMGT/HLA FASTA)
and a batch HMM domain search, used by build_mhc_refdb.py."""
import collections, csv, gzip, re, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from immtools import _mhc as H

csv.field_size_limit(10**8)
CLADES = ("Mammalia", "Aves", "Crocodylia", "Testudines", "Lepidosauria", "Amphibia",
          "Coelacanthiformes", "Dipnoi", "Actinopterygii", "Chondrichthyes")
MIN_SCORE = 40.0


def read_fasta(path):
    name, seq = None, []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = ln[1:].strip(), []
            else:
                seq.append(ln.strip())
    if name:
        yield name, "".join(seq)


def uniprot(seqdir):
    seen = set()
    # UniProt family annotation (MHC class I / II, B2M) + InterPro MHC class II domain
    # annotation (IPR000353 beta N-terminal, IPR001003 alpha N-terminal): many class II beta
    # entries, especially outside mammals, carry only the InterPro annotation
    for f in ("uniprot_mhc1.tsv.gz", "uniprot_mhc2.tsv.gz", "uniprot_b2m.tsv.gz",
              "uniprot_ipr_mhc2b.tsv.gz", "uniprot_ipr_mhc2a.tsv.gz"):
        if not (seqdir / f).exists():
            continue
        for r in csv.DictReader(gzip.open(seqdir / f, "rt"), delimiter="\t"):
            if r["Entry"] in seen:
                continue
            seen.add(r["Entry"])
            lin = r.get("Taxonomic lineage", "")
            clade = next((c for c in CLADES if c in lin), "other")
            yield dict(uid=f"UP:{r['Entry']}", source="UniProt", acc=r["Entry"],
                       species=r["Organism"].split("(")[0].strip(), taxid=r["Organism (ID)"],
                       clade=clade, afdb=bool(r["AlphaFoldDB"].strip()),
                       reviewed=r["Reviewed"] == "reviewed", name=r["Protein names"][:80],
                       seq=r["Sequence"])


def ipd_mhc(seqdir):
    for name, seq in read_fasta(seqdir / "MHC_prot.fasta"):
        m = re.match(r"IPD-MHC:(\S+)\s+(\S+)", name)
        if not m:
            continue
        allele = m.group(2)
        prefix = allele.split("-")[0]
        yield dict(uid=f"IPD:{m.group(1)}", source="IPD-MHC", acc=m.group(1), species=prefix,
                   taxid="", clade="", afdb=False, reviewed=False, name=allele, seq=seq)


def batch_domains(records):
    """hmmsearch every model against all sequences; return per-record hit lists."""
    from pyhmmer import easel, plan7, hmmer
    hmms, meta = H._models()
    abc = hmms[0].alphabet
    clean = [re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", r["seq"].upper()) for r in records]
    block = easel.DigitalSequenceBlock(abc, [
        easel.TextSequence(name=str(i).encode(), sequence=s).digitize(abc)
        for i, s in enumerate(clean)])
    hits = collections.defaultdict(list)
    for hmm, top in zip(hmms, hmmer.hmmsearch(hmms, block, cpus=2, E=1e-3, domE=1e-3)):
        name = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
        for h in top:
            i = int(h.name.decode() if isinstance(h.name, bytes) else h.name)
            for d in h.domains:
                if d.score >= MIN_SCORE:
                    hits[i].append(H._parse_domain(name, d))
    return clean, hits


