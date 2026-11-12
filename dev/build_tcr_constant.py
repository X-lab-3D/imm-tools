"""Build profile HMMs for TCR constant domains (TRAC, TRBC, TRGC, TRDC) with IMGT C-LIKE
labels from IMGT/DomainDisplay (all species, receptor type TR, domain type C), the same
source and parser as the MHC models.

    python dev/build_tcr_constant.py DomainDisplay_TR_C.html immtools/models/tcr_constant
"""
import collections, re, sys
import html as H
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hmmbuild as BH

KIND = {"C-ALPHA": "C_TRAC", "C-BETA-1": "C_TRBC", "C-BETA-2": "C_TRBC", "C-BETA": "C_TRBC",
        "C-GAMMA": "C_TRGC", "C-GAMMA-1": "C_TRGC", "C-GAMMA-2": "C_TRGC", "C-DELTA": "C_TRDC"}


def header_labels(header):
    """BH.header_labels, but an insertion run is based on the standard position right before
    it (IMGT FG loop 111A.., 112A..), not on the preceding tick mark."""
    numline, marker = header
    lab = BH.header_labels(header)
    cols = list(marker)
    i = 0
    while i < len(cols):
        if cols[i] not in ". |" and not (lab[i] or "").startswith("1") or (
                cols[i] not in ". |" and i > 0 and cols[i - 1] in ".|"):
            pass
        i += 1
    # recompute runs
    runs, cur = [], []
    first_bar = cols.index("|")
    for i, c in enumerate(cols):
        if i > first_bar and c not in ". |":
            cur.append(i)
        elif cur:
            runs.append(cur); cur = []
    if cur:
        runs.append(cur)
    for run in runs:
        prev = next((lab[j] for j in range(run[0] - 1, -1, -1)
                     if lab[j] and lab[j].isdigit()), None)
        if prev is None:
            continue
        BH._label_run(lab, [(j, cols[j]) for j in run], int(prev))
    return lab


GRID = None
LET = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def symmetric(res):
    """IMGT rule for the symmetric loops (C-LIKE 84/85 and FG 111/112): k inserted residues
    take 84A.. / 111A.. for the first ceil(k/2) and ..85A / ..112A for the rest, the gap in
    the middle (DomainDisplay sometimes shows the gap at the end of the run)."""
    out = list(res)
    for a in (84, 111):
        idx = [i for i, (l, _) in enumerate(out)
               if l[:-1] in (str(a), str(a + 1)) and l[-1].isalpha()]
        if not idx:
            continue
        k = len(idx); nl = k - k // 2; nr = k // 2
        labs = [f"{a}{LET[j]}" for j in range(nl)] + [f"{a + 1}{LET[j]}" for j in range(nr - 1, -1, -1)]
        for i, l in zip(idx, labs):
            out[i] = (l, out[i][1])
    return out


def merged_grid(label_lists):
    """Column order over all section headers (labels in display order)."""
    lists = [[l for l in ls if l] for ls in label_lists if ls]
    lists.sort(key=len, reverse=True)
    grid = list(dict.fromkeys(lists[0]))
    for ls in lists[1:]:
        for k, l in enumerate(ls):
            if l not in grid:
                prev = next((ls[j] for j in range(k - 1, -1, -1) if ls[j] in grid), None)
                grid.insert(grid.index(prev) + 1 if prev else 0, l)
    return grid


def parse(page):
    """rows dict(kind, species, allele, source, res) with labels from each section's header."""
    events = []
    for m in re.finditer(r"<pre>(.*?)</pre>", page, re.S):
        lines = H.unescape(re.sub(r"<[^>]+>", "", m.group(1))).split("\n")
        for j, ln in enumerate(lines):
            if re.match(r"^\s*[0-9]+\|", ln) and j > 0:
                events.append((m.start(), "hdr", (lines[j - 1], ln)))
                break
    for m in re.finditer(r"<tr>\s*(.*?)</tr>", page, re.S):
        tds = re.findall(r"<td[^>]*>(.*?)</td>", m.group(1), re.S)
        if len(tds) < 9 or "<pre>" not in tds[-1]:
            continue
        sp, grp, allele, acc, dom, dnum, func = [BH._text(x) for x in tds[1:8]]
        if dom == "Domain label":
            continue
        pre = re.search(r"<pre>(.*?)</pre>", tds[-1], re.S).group(1)
        seq = H.unescape(re.sub(r"<[^>]+>", "", pre)).rstrip()
        events.append((m.start(), "rec", dict(species=sp, allele=allele, domain=dom, func=func,
                                              gapped=seq)))
    events.sort(key=lambda e: e[0])
    global GRID
    GRID = merged_grid([header_labels(x) for _, t, x in events if t == "hdr"])
    rows, lab, skipped = [], None, collections.Counter()
    for _, typ, x in events:
        if typ == "hdr":
            try:
                lab = header_labels(x)
            except ValueError:
                lab = None
            continue
        k = KIND.get(x["domain"])
        if k is None or lab is None:
            skipped[x["domain"]] += 1
            continue
        s = x["gapped"]
        if len(s) > len(lab) or any((lab[i] is None) != (s[i] == " ") for i in range(len(s))):
            skipped["grid mismatch"] += 1
            continue
        res = [(l, c) for l, c in zip(lab, s) if l is not None and c not in ". "]
        res = symmetric(res)
        if res:
            rows.append(dict(kind=k, species=x["species"], allele=x["allele"],
                             source="IMGT DomainDisplay", res=res))
    return rows, skipped


if __name__ == "__main__":
    rows, skipped = parse(open(sys.argv[1], errors="ignore").read())
    print(collections.Counter(r["kind"] for r in rows), dict(skipped), file=sys.stderr)
    meta = BH.build_from_rows(rows, {"C": GRID, "G": []}, sys.argv[2])
    meta = {"description": "TCR constant-domain numbering models (IMGT C-LIKE), from "
                           "IMGT/DomainDisplay", "models": {
        k: {kk: v[kk] for kk in ("grid", "n_seqs", "species", "M", "match_labels",
                                 "insert_labels") if kk in v}
        for k, v in meta["models"].items()}}
    import json, pathlib
    (pathlib.Path(sys.argv[2]) / "hmms.json").write_text(json.dumps(meta, indent=1))
