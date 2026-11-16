#!/usr/bin/env python3
"""
TCR-pMHC contacts by IMGT position, pooled over several structures.

Because every structure is IMGT-numbered first, "CDR3 beta 109 touches MHC 1072" means the
same thing in every complex, whatever the species, allele or author numbering. The script
writes a table of contacts and a heat map: TCR positions (CDRs) x MHC positions, coloured by
the number of complexes in which the pair is in contact.

Usage:
    python contact_map.py 1ao7.pdb 1mi5.pdb 3qdg.cif ... [-o contacts]    # -> contacts.tsv, contacts.png

Contact: any protein heavy-atom pair within 4.5 A (change with --cutoff). Give complexes of
one MHC class at a time: positions follow each chain's own IMGT numbering (class I alpha1
1-92, alpha2 1001-1092; see the README).
"""
import argparse
import collections
from pathlib import Path

from Bio.PDB import NeighborSearch
from Bio.PDB.Polypeptide import is_aa

import immtools

CDR = {"CDR1": (27, 38), "CDR2": (56, 65), "CDR3": (105, 117)}
SIDE = {"TCR_alpha": "α", "TCR_beta": "β", "TCR_gamma": "γ", "TCR_delta": "δ"}


def cdr_of(number):
    return next((k for k, (lo, hi) in CDR.items() if lo <= number <= hi), None)


def contacts(path, cutoff=4.5):
    """{(tcr label, mhc label)} for one complex, labels like 'CDR3β 109' and 'MHC 1072'."""
    text = immtools.renumber(path)                        # TCR, MHC, peptide: IMGT numbering
    types = {r["chain"]: r["type"] for r in immtools.annotate(text)}
    cx = immtools.complexes(text)[0]                      # first TCR-pMHC complex
    tcr = [cx[r] for r in ("tcra", "tcrb") if r in cx]
    mhc = [cx[r] for r in ("mhc", "mhcb") if r in cx]
    model = immtools.structure(text)
    atoms = [a for c in tcr + mhc for r in model[c] if is_aa(r)     # protein only, no glycans
             for a in r if a.element != "H"]
    found = set()
    for a, b in NeighborSearch(atoms).search_all(cutoff):
        ra, rb = a.get_parent(), b.get_parent()
        ca, cb = ra.get_parent().id, rb.get_parent().id
        if ca in mhc and cb in tcr:
            ra, rb, ca, cb = rb, ra, cb, ca
        if not (ca in tcr and cb in mhc):
            continue
        num = ra.id[1]
        region = cdr_of(num)
        if region is None:
            continue                                      # framework contacts left out
        t = f"{region}{SIDE.get(types[ca], '?')} {num}{ra.id[2].strip()}"
        m = f"{'MHC' if types[cb] != 'MHC_II_beta' else 'MHCβ'} {rb.id[1]}{rb.id[2].strip()}"
        found.add((t, m))
    return found


def order_tcr(label):
    region, pos = label.split(" ")
    side = "αβγδ".index(region[-1]) if region[-1] in "αβγδ" else 9
    num = int("".join(c for c in pos if c.isdigit()))
    return (side, int(region[3]), num, pos)


def order_mhc(label):
    pos = label.split(" ")[1]
    return (label.startswith("MHCβ"), int("".join(c for c in pos if c.isdigit())), pos)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("structures", nargs="+")
    ap.add_argument("-o", "--out", default="contacts")
    ap.add_argument("--cutoff", type=float, default=4.5)
    a = ap.parse_args()

    count = collections.Counter()
    for f in a.structures:
        pairs = contacts(f, a.cutoff)
        count.update(pairs)
        print(f"{Path(f).name}: {len(pairs)} CDR-MHC contacts")
    with open(a.out + ".tsv", "w") as fh:
        fh.write("tcr_position\tmhc_position\tcomplexes\n")
        for (t, m), n in sorted(count.items(), key=lambda x: (order_tcr(x[0][0]), order_mhc(x[0][1]))):
            fh.write(f"{t}\t{m}\t{n}\n")
    plot(count, len(a.structures), a.out + ".png")
    print(f"wrote {a.out}.tsv and {a.out}.png")


def plot(count, n_structures, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize

    rows = sorted({t for t, _ in count}, key=order_tcr)
    cols = sorted({m for _, m in count}, key=order_mhc)
    grid = [[count.get((t, m), 0) for m in cols] for t in rows]
    # one-hue sequential ramp (blue): 0 = chart surface, darker = more complexes
    ramp = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
    cmap = LinearSegmentedColormap.from_list("blue", ramp)
    norm = Normalize(0, n_structures)

    w = max(6, 0.22 * len(cols) + 2.2)
    h = max(3, 0.22 * len(rows) + 1.6)
    fig, ax = plt.subplots(figsize=(w, h), dpi=200)
    fig.patch.set_facecolor("#fcfcfb")
    ax.imshow(grid, cmap=cmap, norm=norm, aspect="equal")
    ax.set_xticks(range(len(cols)), [c.split(" ")[1] for c in cols], rotation=90, fontsize=6,
                  color="#52514e")
    ax.set_yticks(range(len(rows)), rows, fontsize=6, color="#52514e")
    ax.set_xticks([x - 0.5 for x in range(1, len(cols))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(rows))], minor=True)
    ax.grid(which="minor", color="#fcfcfb", linewidth=1.5)      # 2px surface gap between cells
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlabel("MHC position (IMGT)", fontsize=7, color="#0b0b0b")
    ax.set_ylabel("TCR position (IMGT)", fontsize=7, color="#0b0b0b")
    ax.set_title(f"CDR–MHC contacts in {n_structures} TCR–pMHC complexes", fontsize=8,
                 color="#0b0b0b", loc="left")
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.6,
                      ticks=range(0, n_structures + 1, max(1, n_structures // 5)))
    cb.set_label("complexes with this contact", fontsize=6, color="#52514e")
    cb.ax.tick_params(labelsize=6, colors="#52514e", length=0)
    cb.outline.set_visible(False)
    fig.tight_layout()
    fig.savefig(path, facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.15)


if __name__ == "__main__":
    main()
