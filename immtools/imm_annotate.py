# SPDX-License-Identifier: Apache-2.0
"""
Reports, for each chain, its type, the complex it belongs to, species, closest gene / allele and
domain ranges, as a tab-separated table.

Usage:
    imm_annotate [-<chain id(s)>] <pdb or mmcif file>

Example:
    imm_annotate 1ao7.pdb
    imm_annotate 1ao7.pdb > 1ao7.tsv
    for f in *.pdb; do imm_annotate $f; done | awk 'NR==1 || !/^file/'

Columns:
    file chain type complex length species gene_allele identity domains
    type      TCR_alpha/beta/gamma/delta, MHC_I, MHC_II_alpha, MHC_II_beta, B2M, peptide,
              antibody_heavy/kappa/lambda, other
    complex   TCR-pMHC complex the chain belongs to (1, 2, ... ; see imm_complex)
    TCR       species and V / J germline genes from ANARCI (identity of the V and J region)
    MHC       closest allele / gene by the peptide-binding region (class I alpha1+alpha2,
              class II alpha1 or beta1; beta2m: whole domain) among IPD-IMGT/HLA, IPD-MHC
              and UniProt sequences, with % identity over that region; alleles with the same
              region are listed together ("+N" more)
    domains   name:first-last residue (numbers as in the input file)

This program is part of imm-tools, made in the style of pdb-tools.
"""
import sys

from immtools import _io, _doc, _chains, _complex, _annot, _mhc, _tcr

SPECIES = {"human": "Homo sapiens", "mouse": "Mus musculus", "rat": "Rattus norvegicus",
           "rabbit": "Oryctolagus cuniculus", "rhesus": "Macaca mulatta", "pig": "Sus scrofa",
           "alpaca": "Vicugna pacos", "cow": "Bos taurus"}

COLS = ("file", "chain", "type", "complex", "length", "species", "gene_allele",
        "identity", "domains")


def _range(rs, idx):
    idx = sorted(idx)
    return f"{_io.resid(rs[idx[0]])}-{_io.resid(rs[idx[-1]])}" if idx else ""


def annotate(c):
    rs = c["residues"]
    row = dict(species="", gene_allele="", identity="", domains="")
    if c["mhc"]:
        labels, info = c["mhc"]
        best = _annot.closest(c["seq"], labels, info["chain_type"])
        if best:
            row.update(gene_allele=best[0], species=best[1], identity=f"{best[2]:.1f}")
        parts = []
        for d, m in enumerate(_mhc.CHAIN_TYPES[info["chain_type"]]):
            idx = [i for i, l in enumerate(labels) if l[0] // 1000 == d and 0 < l[0] % 1000 < 200]
            if idx:
                parts.append(f"{m}:{_range(rs, idx)}")
        row["domains"] = ";".join(parts)
    elif c["tcr"]:
        doms = c["tcr"][1]
        d = doms[0]
        row["species"] = SPECIES.get(d["species"], d["species"] or "")
        genes = [g for g in (d.get("v_gene"), d.get("j_gene")) if g]
        row["gene_allele"] = " / ".join(genes)
        ids = [f"{100 * x:.1f}" for x in (d.get("v_identity"), d.get("j_identity")) if x]
        row["identity"] = " / ".join(ids)
        parts = [f"V:{_range(rs, range(x['start'], x['end'] + 1))}" for x in doms]
        lab, info = _tcr.full_labels(c["seq"], constant=True)
        if info["constant"]:
            cb = 1000 * (1 + max(l[0] for l in c["tcr"][0] if l) // 1000)
            ci = [i for i, l in enumerate(lab) if cb < l[0] < cb + 200]
            parts.append(f"C:{_range(rs, ci)}")
        row["domains"] = ";".join(parts)
    return row


def rows(model, chains=None, name=""):
    """One dict per chain (keys: COLS)."""
    ident = _chains.identify_all(model)
    cxs = _complex.find(model, ident)
    where = {cid: k for k, cx in enumerate(cxs, 1) for cid in cx.values()}
    for cid, c in ident.items():
        if c["kind"] == "none" or (chains and cid not in chains):
            continue
        row = dict(file=name, chain=cid, type=c["kind"], length=str(len(c["residues"])),
                   complex=str(where[cid]) if cid in where else "")
        row.update(annotate(c))
        yield row


def run(fhandle, option, name):
    model = _doc.read(fhandle).model()
    yield "\t".join(COLS) + "\n"
    for row in rows(model, _io.chain_list(option), name):
        yield "\t".join(row[k] for k in COLS) + "\n"


def main():
    option, fh = _io.check_input(sys.argv[1:], __doc__)
    _io.chain_list(option)
    _io.write(list(run(fh, option, _io.name_of(fh, "stdin"))))
    fh.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
