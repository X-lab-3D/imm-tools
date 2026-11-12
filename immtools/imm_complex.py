# SPDX-License-Identifier: Apache-2.0
"""
Picks one TCR-pMHC complex out of a file with several copies (crystal asymmetric units).
Chains are grouped by type and contacts: MHC heavy chain (or class II alpha) + beta2m (or
class II beta) + the peptide in its groove + the TCR alpha/beta pair that touches them.
All complexes found are listed on stderr.

Usage:
    imm_complex [-<n>] <pdb or mmcif file>

    n: which complex (1, 2, ...; default 1), in the order of their first chain in the file.

Example:
    imm_complex 3qdg.pdb > complex1.pdb
    imm_complex -2 3qdg.pdb | imm_rechain > complex2_standard.pdb

TCRs without pMHC, and lone MHC chains, form complexes of their own.

This program is part of imm-tools, made in the style of pdb-tools.
"""
import sys

from immtools import _io, _doc, _complex


def run(fhandle, option):
    n = int(option) if option else 1
    doc = _doc.read(fhandle)
    model = doc.model()
    cxs = _complex.find(model)
    for k, cx in enumerate(cxs, 1):
        sys.stderr.write(f"complex {k}: {_complex.describe(cx)}\n")
    if not cxs:
        sys.stderr.write("ERROR!! No TCR, MHC or beta2m chain found.\n")
        sys.exit(1)
    if not 1 <= n <= len(cxs):
        sys.stderr.write(f"ERROR!! Complex {n} requested, {len(cxs)} found.\n")
        sys.exit(1)
    return doc.selected(set(cxs[n - 1].values())).text()


def main():
    option, fh = _io.check_input(sys.argv[1:], __doc__)
    if option and not option.isdigit():
        sys.stderr.write(f"ERROR!! Option must be a complex number: '-{option}'\n")
        sys.stderr.write(__doc__)
        sys.exit(1)
    _io.write(run(fh, option))
    fh.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
