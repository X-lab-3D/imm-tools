#!/usr/bin/env python3
"""
CDR3 table of many TCR structures: one row per TCR chain with V / J genes, CDR3 sequence and
its residues in IMGT numbering, so loops of different length line up (IMGT inserts the extra
residues of a long CDR3 symmetrically at 111A, 112A, ...).

Usage:
    python cdr3_table.py 1ao7.pdb 1mi5.pdb 3qdg.cif ... > cdr3.tsv
"""
import sys
from pathlib import Path

import immtools


def main(paths):
    print("file\tchain\ttype\tv_gene\tj_gene\tcdr3\tlength\timgt_positions")
    for path in paths:
        text = immtools.renumber(path, mhc=False, peptide=False)   # TCR chains only
        genes = {r["chain"]: r["gene_allele"].split(" / ") for r in immtools.annotate(text)}
        for c in immtools.cdrs(text):
            if c["cdr"] != "CDR3":
                continue
            v, j = (genes.get(c["chain"], []) + ["", ""])[:2]
            positions = ",".join(r.strip() for r in c["residues"])
            print(f"{Path(path).stem}\t{c['chain']}\t{c['type']}\t{v}\t{j}\t{c['sequence']}\t"
                  f"{c['length']}\t{positions}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
