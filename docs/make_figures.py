#!/usr/bin/env python3
"""Make the README figure: the alpha1 helix of MHC class I from five species, IMGT-numbered.

    python docs/make_figures.py          # -> docs/images/alpha1_species.png

The same IMGT position is the same place in the structure in every species, although the
sequence positions differ (chicken and carp have an insertion, 54A). The bars on top show in
how many of 10 human and mouse TCR-pMHC complexes the position touches the TCR
(examples/contact_map.py on 1AO7 1MI5 3QDG 1OGA 5BRZ 2VLJ 3GSN 4FTV 3HG1 2CKB).
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import immtools

# MHC class I heavy chains (alpha1 + alpha2) of PDB entries
CHAINS = [
    ("human HLA-A*02:01", "1AO7 A",
     "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEA"
     "GSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAP"),
    ("mouse H-2Kb", "1G6R H",
     "GPHSLRYFVTAVSRPGLGEPRYMEVGYVDDTEFVRFDSDAENPRYEPRARWMEQEGPEYWERETQKAKGNEQSFRVDLRTLLGYYNQSKG"
     "GSHTIQVISGCEVGSDGRLLRGYQQYAYDGCDYIALNEDLKTWTAADMAALITKHKWEQAGEAERLRAYLEGTCVEWLRRYLKNGNATLLRTDSP"),
    ("pig SLA-1*04:01", "3QQ3 A",
     "GPHSLSYFYTAVSRPDRGDSRFIAVGYVDDTQFVRFDNYAPNPRMEPRVPWIQQEGQEYWDRETRNVKETAQTYGVGLNTLRGYYNQSEA"
     "GSHTLQSMYGCYLGPDGLLLHGYRQDAYDGADYIALNEDLRSWTAADMAAQITKRKWEAADEAERRRSYLQGLCVESLRRYLEMGKDTLQRAEPP"),
    ("chicken BF2*21:01", "3BEV A",
     "ELHTLRYIRTAMTDPGPGLPWFVDVGYVDGELFMHYNSTARRAVPRTEWIAANTDQQYWDRETQIVQGSEQINRENLDILRRRYNQTG"
     "GSHTVQWMSGCDILEDGTIRGYHQAAYDGRDFVAFDKGTMTLTAAVPEAVPTKRKWEEGGYAEGLKQYLEETCVEWLRRYVEYGKAELGRRERPEVR"),
    ("grass carp UAA", "5CNZ A",
     "GTHSLKYVYTGVSRGIDFPEFTAVGMVDDGQFMYFDSNSMKAVPKTEWIRQNEGADYWDRQTQVLIGAHQVFKDSIQIVMERFNQSKG"
     "VHTWQNMYGCELNDDGTTQGFYQYAYDGEDFVSLDKNTLTWTAANPQAVITKHKWEALAVAEQNKGYLENTCIEWLKKYVAYGKDTLERKVSPQVSL"),
]
# complexes (of 10) in which the alpha1 position contacts the TCR (examples/contact_map.py)
TCR_CONTACTS = {55: 1, 58: 4, 59: 1, 62: 6, 65: 8, 66: 10, 68: 5, 69: 8, 70: 2, 72: 8, 73: 9,
                75: 3, 76: 7, 79: 2, 80: 1}
FIRST, LAST = 50, 92

INK, INK2, MUTED, SURFACE = "#0b0b0b", "#52514e", "#898781", "#fcfcfb"
BLUE, TINT, INSERT = "#2a78d6", "#e7f0fc", "#fde7dd"


def main(out=Path(__file__).parent / "images" / "alpha1_species.png"):
    rows = []
    for name, pdb, seq in CHAINS:
        num = immtools.number_mhc(seq)
        pos = {}
        for i, ((n, ic), aa) in enumerate(zip(num.numbers, seq)):
            if FIRST <= n <= LAST:
                pos[(n, ic)] = (aa, i + 1)          # residue, sequence position
        rows.append((name, pdb, pos))
    cols = sorted({k for _, _, p in rows for k in p}, key=lambda k: (k[0], k[1]))

    w, h = len(cols), len(rows)
    fig = plt.figure(figsize=(0.2 * w + 2.6, 2.75), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    bars = fig.add_axes([0.17, 0.66, 0.70, 0.17])
    grid = fig.add_axes([0.17, 0.12, 0.70, 0.52], sharex=bars)
    for ax in (bars, grid):
        ax.set_facecolor(SURFACE)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)

    # TCR contacts per position
    x = range(w)
    v = [TCR_CONTACTS.get(n, 0) if not ic else 0 for n, ic in cols]
    bars.bar(x, v, width=0.6, color=BLUE)
    bars.axhline(0, color="#c3c2b7", linewidth=0.6)
    bars.set_ylim(0, 10)
    bars.set_yticks([0, 5, 10])
    bars.tick_params(labelsize=5, colors=MUTED, labelbottom=False)
    bars.set_title("TCR contacts (of 10 complexes)", fontsize=6, color=INK2, loc="left", pad=3)
    for xi, val, (n, _) in zip(x, v, cols):
        if val >= 8:
            bars.text(xi, val + 0.4, str(n), ha="center", va="bottom", fontsize=4.5, color=INK)

    # residues
    for xi, (n, ic) in enumerate(cols):
        if ic:
            grid.add_patch(Rectangle((xi - 0.5, -0.5), 1, h, color=INSERT, lw=0))
        elif TCR_CONTACTS.get(n, 0) >= 5:
            grid.add_patch(Rectangle((xi - 0.5, -0.5), 1, h, color=TINT, lw=0))
    for yi, (name, pdb, pos) in enumerate(rows):
        for xi, c in enumerate(cols):
            aa = pos.get(c, ("–", None))[0]
            grid.text(xi, yi, aa, ha="center", va="center", fontsize=6.5, family="monospace",
                      color=INK if aa != "–" else MUTED)
        seq65 = pos.get((65, ""), (None, "?"))[1]
        grid.text(w - 0.2, yi, f"{seq65}", ha="left", va="center", fontsize=5.5, color=INK2)
    grid.text(w - 0.2, -1.0, "seq. no.\nof 65", ha="left", va="bottom", fontsize=4.5,
              color=MUTED)
    grid.set_xlim(-0.5, w - 0.5)
    grid.set_ylim(h - 0.5, -0.5)
    grid.set_yticks(range(h), [f"{n}  ({p})" for n, p, _ in rows], fontsize=5.5, color=INK2)
    grid.set_xticks(range(w), [f"{n}{ic}" if (n % 5 == 0 or ic) else "" for n, ic in cols],
                    fontsize=5, color=INK2)
    grid.set_xlabel("IMGT position (G-ALPHA1 domain, helix)", fontsize=6, color=INK, labelpad=3)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, facecolor=SURFACE, bbox_inches="tight", pad_inches=0.08)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
