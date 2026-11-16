#!/usr/bin/env python3
"""
Write a PyMOL script that colours a TCR-pMHC structure by IMGT region: CDR1/2/3 of each TCR
chain, the two helices of the MHC groove, the peptide. Because it works on IMGT numbers, the
same script colours any complex the same way.

Usage:
    python pymol_regions.py [-footprint] 3qdg_imgt.pdb > 3qdg_regions.pml
    pymol 3qdg_imgt.pdb 3qdg_regions.pml

-footprint: show only the peptide-binding groove, the peptide and the CDR loops, looking down
from the TCR onto the MHC (the TCR "footprint").

The input must be IMGT-numbered (imm_tcrnum | imm_mhcnum | imm_pepnum, or immtools.renumber).
Regions (IMGT): CDR1 27-38, CDR2 56-65, CDR3 105-117; groove helices 50-92 (class I alpha1,
class II alpha) and 1050-1092 (class I alpha2, class II beta: 50-92 of its G domain).
"""
import sys
from pathlib import Path

import immtools

# categorical palette, fixed order; text and frameworks stay neutral
COLOURS = {"cdr1": "0xeda100", "cdr2": "0x1baf7a", "cdr3": "0xeb6834",
           "helix": "0x2a78d6", "peptide": "0xe87ba4",
           "tcr": "0xc3c2b7", "mhc": "0xe1e0d9", "other": "0xf0efec"}
CDR = {"cdr1": (27, 38), "cdr2": (56, 65), "cdr3": (105, 117)}
HELIX = {"MHC_I": [(50, 92), (1050, 1092)], "MHC_II_alpha": [(50, 92)],
         "MHC_II_beta": [(50, 92)]}


def residues(chain, ranges):
    """PyMOL residue list (with insertion codes) of the residues numbered within `ranges`."""
    out = [f"{r.id[1]}{r.id[2].strip()}" for r in chain
           if r.id[0] == " " and any(lo <= r.id[1] <= hi for lo, hi in ranges)]
    return "+".join(out)


FOOTPRINT = """\
hide cartoon, {obj} and not ({groove} or {cdrs} or {pep})
orient {obj} and ({groove})
python
from pymol import cmd
g = cmd.centerofmass("{obj} and ({groove})")
t = cmd.centerofmass("{obj} and ({cdrs})")
v = cmd.get_view()
if sum((t[i] - g[i]) * v[3 * i + 2] for i in range(3)) < 0:   # TCR behind: turn round
    cmd.turn("x", 180)
python end
"""


def main(path, footprint=False):
    obj = Path(path).stem
    model = immtools.structure(path)
    types = {r["chain"]: r["type"] for r in immtools.annotate(path)}
    lines = [f"set_color imm_{k}, [{int(v[2:4], 16)}, {int(v[4:6], 16)}, {int(v[6:8], 16)}]"
             for k, v in COLOURS.items()]
    lines += ["hide everything", f"show cartoon, {obj}", f"color imm_other, {obj}"]

    groove, cdrs, pep = [], [], []

    def colour(name, chain, ranges, c):
        res = residues(chain, ranges)
        if res:
            lines.append(f"select {name}, {obj} and chain {chain.id} and resi {res}")
            lines.append(f"color imm_{c}, {name}")
            if c.startswith("cdr"):
                cdrs.append(name)

    for ch in model:
        t = types.get(ch.id, "")
        if t.startswith("TCR_"):
            side = t[4]                                      # a, b, g, d
            lines.append(f"color imm_tcr, {obj} and chain {ch.id}")
            for k, rng in CDR.items():
                colour(f"{k}{side}", ch, [rng], k)
        elif t in HELIX:
            lines.append(f"color imm_mhc, {obj} and chain {ch.id}")
            colour(f"helix_{ch.id}", ch, HELIX[t], "helix")
            g = "1-1092" if t == "MHC_I" else "1-92"           # G domain(s)
            groove.append(f"(chain {ch.id} and resi {g})")
        elif t == "peptide":
            pep.append(f"chain {ch.id}")
            lines += [f"color imm_peptide, {obj} and chain {ch.id}",
                      f"show sticks, {obj} and chain {ch.id} and not name N+C+O"]
    if footprint and groove and cdrs:
        lines += FOOTPRINT.format(obj=obj, groove=" or ".join(groove), cdrs=" or ".join(cdrs),
                                  pep=" or ".join(pep) or "none").splitlines()
    else:
        lines.append("orient")
    lines += ["deselect", ""]
    sys.stdout.write("\n".join(lines))


if __name__ == "__main__":
    args = sys.argv[1:]
    fp = "-footprint" in args
    args = [a for a in args if a != "-footprint"]
    if len(args) != 1:
        sys.exit(__doc__)
    main(args[0], fp)
