# SPDX-License-Identifier: Apache-2.0
"""
Writes the TCR, the pMHC, the MHC (without peptide) and the peptide of each TCR-pMHC
complex in the file to separate files (same format as the input):

    <name>_TCR.pdb  <name>_pMHC.pdb  <name>_MHC.pdb  <name>_peptide.pdb   (.cif for mmCIF)

(<name>_<n>_... when the file holds several complexes; <name> is the input file name, or
'split' for stdin). MHC includes beta2m. The files are written to the current directory
and listed on stderr; nothing goes to stdout.

Usage:
    imm_split [-<parts>] <pdb or mmcif file>

    parts: comma-separated subset of tcr, pmhc, mhc, pep (default: all four)

Example:
    imm_split 1ao7.pdb
    imm_split -tcr,pmhc 1ao7.pdb
    imm_complex -1 3qdg.pdb | imm_split

This program is part of imm-tools, made in the style of pdb-tools.
"""
import sys

from immtools import _io, _doc, _complex
from immtools._registry import PARTS


def parts_of(option):
    parts = [p.strip().lower() for p in option.split(",")] if option else list(PARTS)
    bad = [p for p in parts if p not in PARTS]
    if bad:
        sys.stderr.write(f"ERROR!! Unknown part(s): {', '.join(bad)}\n")
        sys.stderr.write(__doc__)
        sys.exit(1)
    return parts


def run(fhandle, option, name):
    parts = parts_of(option)
    doc = _doc.read(fhandle)
    model = doc.model()
    cxs = [cx for cx in _complex.find(model) if len(cx) > 1]
    if not cxs:
        sys.stderr.write("ERROR!! No complex found.\n")
        sys.exit(1)
    written = []
    for k, cx in enumerate(cxs, 1):
        stem = name if len(cxs) == 1 else f"{name}_{k}"
        for p in parts:
            label, roles = PARTS[p]
            keep = {cx[r] for r in roles if r in cx}
            if not keep:
                continue
            out = f"{stem}_{label}{doc.ext}"
            with open(out, "w") as fh:
                fh.write(doc.selected(keep).text())
            written.append(out)
            sys.stderr.write(f"{out}: chains {','.join(sorted(keep))}\n")
    return written


def main():
    option, fh = _io.check_input(sys.argv[1:], __doc__)
    parts_of(option)
    run(fh, option, _io.name_of(fh, "split"))
    fh.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
