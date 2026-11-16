# SPDX-License-Identifier: Apache-2.0
"""
Renames the chains of a TCR-pMHC complex to a standard layout:

    A  TCR alpha (or gamma)        M  MHC class I heavy chain / class II alpha
    B  TCR beta (or delta)         N  MHC class II beta
    P  peptide                     2  beta2-microglobulin

Other chains keep their identifier unless it is one of the new ones (then the next free
letter, reported on stderr). The file must hold one complex: use imm_complex first for
files with several copies.

Usage:
    imm_rechain [-<layout>] <pdb or mmcif file>

    layout: role:ID pairs for the chains you want named differently; roles tcra, tcrb,
            mhc, mhcb (class II beta), b2m, pep. E.g. -mhc:A,b2m:B,pep:C,tcra:D,tcrb:E

Example:
    imm_rechain 1ao7.pdb > 1ao7_std.pdb
    imm_rechain -mhc:H,b2m:L 1ao7.pdb            # other letters for MHC and beta2m
    imm_complex -1 3qdg.pdb | imm_rechain | imm_tcrnum | imm_mhcnum | imm_pepnum

This program is part of imm-tools, made in the style of pdb-tools.
"""
import string
import sys

from immtools import _io, _doc, _complex
from immtools._registry import LAYOUT as DEFAULT, ROLES


def layout(option):
    lay = dict(DEFAULT)
    if not option:
        return lay
    parts = [p.strip() for p in option.split(",")]
    if all(":" in p for p in parts):
        for p in parts:
            r, c = p.split(":", 1)
            if r.lower() not in lay or len(c) != 1:
                sys.stderr.write(f"ERROR!! Invalid layout entry: '{p}'\n")
                sys.exit(1)
            lay[r.lower()] = c
    else:
        sys.stderr.write(f"ERROR!! Invalid layout: '-{option}'\n")
        sys.stderr.write(__doc__)
        sys.exit(1)
    return lay


def run(fhandle, option):
    lay = layout(option)
    doc = _doc.read(fhandle)
    model = doc.model()
    cxs = [cx for cx in _complex.find(model) if len(cx) > 1] or _complex.find(model)
    if len(cxs) != 1:
        sys.stderr.write(f"ERROR!! {len(cxs)} complexes in the file; select one with "
                         "imm_complex -<n> first.\n")
        for k, cx in enumerate(cxs, 1):
            sys.stderr.write(f"  complex {k}: {_complex.describe(cx)}\n")
        sys.exit(1)
    cx = cxs[0]
    mapping = {cid: lay[r] for r, cid in cx.items()}
    targets = set(mapping.values())
    if len(targets) != len(mapping):
        sys.stderr.write("ERROR!! Two roles of this complex map to the same chain ID.\n")
        sys.exit(1)
    others = [c.id for c in model if c.id not in mapping]
    taken = targets | {c for c in others if c not in targets}
    free = [c for c in string.ascii_uppercase + string.digits + string.ascii_lowercase
            if c not in taken]
    for c in others:
        if c in targets:
            mapping[c] = free.pop(0)
            sys.stderr.write(f"chain {c} (not in the complex) renamed to {mapping[c]}\n")
    for r in ROLES:
        if r in cx:
            sys.stderr.write(f"{r}: {cx[r]} -> {lay[r]}\n")
    return doc.renamed(mapping).text()


def main():
    option, fh = _io.check_input(sys.argv[1:], __doc__)
    layout(option)
    _io.write(run(fh, option))
    fh.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
