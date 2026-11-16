# SPDX-License-Identifier: Apache-2.0
"""
Reports, for each chain, its type, the complex it belongs to, species, closest gene / allele and
domain ranges, as a tab-separated table.

Usage:
    imm_annotate [-<chain id(s)>,freq:<source>[/<population>]] <pdb or mmcif file>

Example:
    imm_annotate 1ao7.pdb
    imm_annotate 1ao7.pdb > 1ao7.tsv
    imm_annotate -D,E 1ao7.pdb                   # chains D and E only
    imm_annotate -freq:vdjonline 1ao7.pdb        # allele-frequency source (default vdjonline)
    for f in *.pdb; do imm_annotate $f; done | awk 'NR==1 || !/^file/'

Columns:
    file chain type complex length species gene_allele identity closest_common_allele
    common_identity common_differences common_source domains
    type      TCR_alpha/beta/gamma/delta, MHC_I, MHC_II_alpha, MHC_II_beta, B2M, peptide,
              antibody_heavy/kappa/lambda, other
    complex   TCR-pMHC complex the chain belongs to (1, 2, ... ; see imm_complex)
    TCR       species and V / J germline genes from ANARCI (identity of the V and J region)
    closest_common_allele, common_*
              TCRs: the closest common allele of the V and J gene (frequency >= 5% in the
              allele-frequency source), its identity over the germline-encoded positions
              (V: IMGT 1-104, J: 118-128), the differences from it, e.g. S57N (germline Ser
              at IMGT 57, Asn in this chain; "none" if identical) and the source used.
              Sources and populations: README, "Annotation". Default: vdjonline (human; TCR
              repertoires of 134 donors, predominantly of European descent)
    domains   name:first-last residue (numbers as in the input file)

This program is part of imm-tools, made in the style of pdb-tools.
"""
import sys

from immtools import _io, _doc, _chains, _complex, _annot, _mhc, _tcr, _alleles

SPECIES = {"human": "Homo sapiens", "mouse": "Mus musculus", "rat": "Rattus norvegicus",
           "rabbit": "Oryctolagus cuniculus", "rhesus": "Macaca mulatta", "pig": "Sus scrofa",
           "alpaca": "Vicugna pacos", "cow": "Bos taurus"}

COLS = ("file", "chain", "type", "complex", "length", "species", "gene_allele",
        "identity", "closest_common_allele", "common_identity", "common_differences",
        "common_source", "domains")


def _range(rs, idx):
    idx = sorted(idx)
    return f"{_io.resid(rs[idx[0]])}-{_io.resid(rs[idx[-1]])}" if idx else ""


def options(option):
    """'-D,E,freq:pmtr/EUR' -> ({'D', 'E'}, 'pmtr/EUR'); checks the source."""
    parts = [p.strip() for p in option.split(",") if p.strip()]
    freq = [p.split(":", 1)[1] for p in parts if p.startswith("freq:")]
    chains = _io.chain_list(",".join(p for p in parts if not p.startswith("freq:")))
    choice = freq[-1] if freq else ""
    try:
        _alleles.check(choice)
    except ValueError as e:
        sys.stderr.write(f"ERROR!! {e}\n")
        sys.exit(1)
    return chains, choice


def annotate(c, freq=""):
    rs = c["residues"]
    row = dict(species="", gene_allele="", identity="", domains="", **_alleles.EMPTY)
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
        row.update(_alleles.annotate(c["seq"], c["tcr"][0], d, freq))
        parts = [f"V:{_range(rs, range(x['start'], x['end'] + 1))}" for x in doms]
        lab, info = _tcr.full_labels(c["seq"], constant=True)
        if info["constant"]:
            cb = 1000 * (1 + max(l[0] for l in c["tcr"][0] if l) // 1000)
            ci = [i for i, l in enumerate(lab) if cb < l[0] < cb + 200]
            parts.append(f"C:{_range(rs, ci)}")
        row["domains"] = ";".join(parts)
    return row


def rows(model, chains=None, name="", freq=""):
    """One dict per chain (keys: COLS). freq: allele-frequency source[/population]."""
    ident = _chains.identify_all(model)
    cxs = _complex.find(model, ident)
    where = {cid: k for k, cx in enumerate(cxs, 1) for cid in cx.values()}
    for cid, c in ident.items():
        if c["kind"] == "none" or (chains and cid not in chains):
            continue
        row = dict(file=name, chain=cid, type=c["kind"], length=str(len(c["residues"])),
                   complex=str(where[cid]) if cid in where else "")
        row.update(annotate(c, freq))
        yield row


def run(fhandle, option, name):
    chains, freq = options(option)
    model = _doc.read(fhandle).model()
    yield "\t".join(COLS) + "\n"
    for row in rows(model, chains, name, freq):
        yield "\t".join(row[k] for k in COLS) + "\n"


def main():
    option, fh = _io.check_input(sys.argv[1:], __doc__)
    options(option)                   # fail early on a bad option
    _io.write(list(run(fh, option, _io.name_of(fh, "stdin"))))
    fh.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
