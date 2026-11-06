# SPDX-License-Identifier: Apache-2.0
"""Helpers for the model-building scripts in dev/: IMGT/DomainDisplay header parsing and
profile-HMM construction from IMGT-labelled rows (one HMM per domain kind; match columns
are the grid columns occupied in at least MATCH_OCC of the sequences)."""
import collections, json, re, sys
import html as H
from pathlib import Path

MATCH_OCC = 0.5
LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def _text(s):
    return H.unescape(re.sub(r"<[^>]+>", "", s)).strip()


def _ticks(numline, n_bars):
    """Tick values in order; adjacent ticks may be printed glued together ('2326')."""
    vals = []
    for tok in re.findall(r"\d+", numline):
        if len(tok) == 4 and len(vals) + 2 <= n_bars and int(tok[:2]) < int(tok[2:]):
            vals += [int(tok[:2]), int(tok[2:])]
        elif len(tok) == 6:
            vals += [int(tok[:3]), int(tok[3:])]
        else:
            vals.append(int(tok))
    if len(vals) != n_bars:
        raise ValueError(f"{len(vals)} tick numbers for {n_bars} '|' markers")
    return vals


def header_labels(header, kind=None):
    """Column index -> IMGT label ('7A', '49', '1C', ...) or None for separator columns.
    '|' columns carry the tick number above them; '.' columns between two ticks are the
    standard positions in between (if there are fewer columns than positions, IMGT leaves
    the central positions out, as in the C-LIKE BC loop 27-38); digits/letters are
    insertion positions."""
    numline, marker = header
    cols = list(marker)
    labels = [None] * len(cols)
    bars = [i for i, c in enumerate(cols) if c == "|"]
    ticks = dict(zip(bars, _ticks(numline, len(bars))))
    # N-terminal insertion digits before position 1: '0987654321' -> 1J..1A
    for i in range(bars[0]):
        if cols[i].isdigit():
            labels[i] = f"1{LETTERS[(int(cols[i]) or 10) - 1]}"
    for b0, b1 in zip(bars, bars[1:] + [None]):
        p0 = ticks[b0]
        labels[b0] = str(p0)
        seg = list(range(b0 + 1, b1 if b1 is not None else len(cols)))
        dots = [i for i in seg if cols[i] == "."]
        ins = [i for i in seg if cols[i] not in ". |"]
        if b1 is None:
            std = list(range(p0 + 1, p0 + 1 + len(dots)))
        else:
            std = list(range(p0 + 1, ticks[b1]))
            extra = len(std) - len(dots)
            if extra > 0:                       # drop central positions
                mid = len(std) // 2
                lo = mid - extra // 2
                std = std[:lo] + std[lo + extra:]
            elif extra < 0:
                raise ValueError(f"more columns than positions between {p0} and {ticks[b1]}")
        for i, p in zip(dots, std):
            labels[i] = str(p)
        if ins:
            _label_run(labels, [(i, cols[i]) for i in ins], p0)
    return labels


def _label_run(labels, run, prev_pos):
    """Insertion run after standard position prev_pos.  Digits '1234567' -> prevA..; letters
    'ABC' -> prevA..; a palindromic digit run '12345677654321' (C-LIKE 84/85) -> 84A..84G
    then 85G..85A."""
    chars = "".join(c for _, c in run)
    if chars.isdigit():
        vals = [int(c) for c in chars]
        peak = vals.index(max(vals))
        split = peak + 1
        if split < len(vals) and vals[split] == vals[peak]:      # '...77...' symmetric
            for k, (i, c) in enumerate(run):
                if k <= peak:
                    labels[i] = f"{prev_pos}{LETTERS[int(c) - 1]}"
                else:
                    labels[i] = f"{prev_pos + 1}{LETTERS[int(c) - 1]}"
            return
        for i, c in run:
            labels[i] = f"{prev_pos}{LETTERS[int(c) - 1]}"
    else:
        for i, c in run:
            labels[i] = f"{prev_pos}{c}"


# --------------------------------------------------------------------------- #
#  assemble per-kind alignments
# --------------------------------------------------------------------------- #
def build_from_rows(rows, grid, outdir, match_occ=MATCH_OCC):
    from pyhmmer import easel, plan7
    abc = easel.Alphabet.amino()
    groups = collections.defaultdict(list)
    for r in rows:
        col_labels = grid["G" if r["kind"].startswith("G_") else "C"]
        idx = {l: j for j, l in enumerate(col_labels)}
        cols = ["-"] * len(col_labels)
        ok = False
        for l, c in r["res"]:
            if l in idx:
                cols[idx[l]] = c
                ok = True
        if ok:
            groups[r["kind"]].append((r, "".join(cols)))

    meta = {"source": "IMGT/3Dstructure-DB DomainDisplay (https://www.imgt.org) + "
                      "structure-checked IMGT-numbered structures",
            "match_occupancy": match_occ, "models": {}}
    hmms = []
    builder = plan7.Builder(abc, architecture="hand")
    bg = plan7.Background(abc)
    for k, items in sorted(groups.items()):
        g = "G" if k.startswith("G_") else "C"
        col_labels = grid[g]
        seen, seqs, species, sources = set(), [], collections.Counter(), collections.Counter()
        for r, cols in items:
            if cols in seen or set(cols) <= {"-"}:
                continue
            seen.add(cols)
            seqs.append(cols)
            species[r["species"]] += 1
            sources[r["source"]] += 1
        n = len(seqs)
        # occupancy over the sequences that actually cover a column (partial alleles are
        # common in IMGT: leading/trailing gaps mean 'unknown', not 'deleted')
        span = []
        for sq in seqs:
            res = [j for j, c in enumerate(sq) if c != "-"]
            span.append((res[0], res[-1]))
        occ = []
        for j in range(len(col_labels)):
            cov = [sq for sq, (a, b) in zip(seqs, span) if a <= j <= b]
            occ.append(sum(sq[j] != "-" for sq in cov) / len(cov) if cov else 0.0)
        is_match = [o >= match_occ for o in occ]
        msa = easel.TextMSA(name=k.encode(), sequences=[
            easel.TextSequence(name=f"{k}_{i}".encode(), sequence=s) for i, s in enumerate(seqs)])
        msa.reference = "".join("x" if m else "." for m in is_match)
        hmm, _, _ = builder.build_msa(msa.digitize(abc), bg)
        hmm.name = k.encode()
        hmms.append(hmm)
        match_labels = [l for l, m in zip(col_labels, is_match) if m]
        between = collections.defaultdict(list)       # key: index of preceding match (0 = N-term)
        m_idx = 0
        for l, m in zip(col_labels, is_match):
            if m:
                m_idx += 1
            else:
                between[m_idx].append(l)
        meta["models"][k] = {
            "grid": g, "n_seqs": n, "species": dict(species), "sources": dict(sources),
            "M": hmm.M, "match_labels": match_labels,
            "insert_labels": {str(i): v for i, v in between.items()},
        }
        print(f"  {k:9s} {n:5d} unique seqs, {len(species):2d} species, M={hmm.M}",
              file=sys.stderr)
        assert hmm.M == len(match_labels)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "mhc_hmms.hmm", "wb") as fh:
        for h in hmms:
            h.write(fh)
    (outdir / "mhc_hmms.json").write_text(json.dumps(meta, indent=1))
    return meta


