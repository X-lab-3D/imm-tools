# SPDX-License-Identifier: Apache-2.0
"""
Prints the CDRs of the TCR chains (IMGT definition: CDR1 27-38, CDR2 56-65, CDR3 105-117)
as a table, or selects CDR residues from the structure.

Usage:
    imm_cdr [-<cdr(s)>] <pdb or mmcif file>

    no option     table: file, chain, type, cdr, residues (input numbering), length, sequence
    cdr(s)        coordinates of the selected CDRs only (PDB or mmCIF, as the input):
                  cdr1, cdr2, cdr3 (both chains), a1, a2, a3 (alpha / gamma chain),
                  b1, b2, b3 (beta / delta chain)

Example:
    imm_cdr 1ao7.pdb
    imm_cdr -cdr3 1ao7.pdb > cdr3.pdb
    imm_cdr -a3,b3 1ao7.pdb | pdb_tofasta

The CDRs are found by IMGT-numbering the chains internally (ANARCI), so the input can have
any residue numbering; the output keeps the input numbering.

This program is part of imm-tools, made in the style of pdb-tools.
"""
import sys

from immtools import _io, _doc, _chains, _tcr
from immtools._registry import CDR, CDR_SIDE as SIDE


def selection(option):
    sel = set()
    for t in [x.strip().lower() for x in option.split(",") if x.strip()]:
        if t in ("cdr1", "cdr2", "cdr3"):
            sel |= {(s, int(t[3])) for s in set(SIDE.values())}
        elif len(t) == 2 and t[0] in set(SIDE.values()) and t[1] in "123":
            sel.add((t[0], int(t[1])))
        else:
            sys.stderr.write(f"ERROR!! Unknown CDR: '{t}'\n")
            sys.stderr.write(__doc__)
            sys.exit(1)
    return sel


def cdrs(model):
    """[(chain, kind, cdr number, [residues])]"""
    out = []
    for ch in model:
        c = _chains.identify(ch)
        if c["kind"] not in SIDE:
            continue
        labels = c["tcr"][0]
        for k, (lo, hi) in CDR.items():
            rs = [r for r, l in zip(c["residues"], labels) if l and l[0] < 1000 and lo <= l[0] <= hi]
            out.append((ch.id, c["kind"], k, rs))
    return out


def run(fhandle, option, name):
    doc = _doc.read(fhandle)
    model = doc.model()
    found = cdrs(model)
    if not option:
        from immtools._chains import sequence
        yield "file\tchain\ttype\tcdr\tresidues\tlength\tsequence\n"
        for cid, kind, k, rs in found:
            rng = f"{_io.resid(rs[0])}-{_io.resid(rs[-1])}" if rs else ""
            yield f"{name}\t{cid}\t{kind}\tCDR{k}\t{rng}\t{len(rs)}\t{sequence(rs)}\n"
        return
    sel = selection(option)
    keep = {_io.residue_key(r) for cid, kind, k, rs in found if (SIDE[kind], k) in sel for r in rs}
    yield doc.residues(keep).text()

def main():
    option, fh = _io.check_input(sys.argv[1:], __doc__)
    if option:
        selection(option)
    _tcr.require()
    _io.write(list(run(fh, option, _io.name_of(fh, "stdin"))))
    fh.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
