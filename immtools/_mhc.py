#!/usr/bin/env python3
"""
_mhc.py (imm-tools) - IMGT unique numbering of MHC chains with profile HMMs (pyhmmer).

Models (models/mhc_broad and models/mhc_mammal, one HMM per domain kind):
    G_ALPHA1, G_ALPHA2, C_ALPHA3   class I heavy chain (alpha1, alpha2, alpha3)
    C_B2M                          beta2-microglobulin
    G_ALPHA, C_IIA                 class II alpha (alpha1, alpha2)
    G_BETA,  C_IIB                 class II beta  (beta1, beta2)

Each HMM match state carries an IMGT label (G-DOMAIN grid 1-92 with 7A, 49A-G, 54A-C, 61A/B,
72A, N-terminal 1J..1A; C-LIKE grid 1-128 with 15A-C, 45A-G, 84A-G/85G-A, 96A/B, 1H..1A).
Residues aligned to match states take that label; residues in insert states take the IMGT
insertion positions defined between those match states (symmetric fill for the 84/85 loop),
and if there are more residues than defined positions, extra letters are appended.

Numbering of a chain:  domain i of the chain type gets offset 1000*i
    class I heavy: alpha1 +0, alpha2 +1000, alpha3 +2000      beta2m: +0
    class II alpha: alpha1 +0, alpha2 +1000                    class II beta: beta1 +0, beta2 +1000
    residues before the first domain: 0, -1, ...;   after the last domain: next free thousand, +1, ...
    residues between domains (linkers) continue the preceding domain: last label + letters
"""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODEL_DIR = None                   # set to one model directory to use that set alone
LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

CHAIN_TYPES = {                     # chain type -> ordered domain models
    "MHC_I_ALPHA": ["G_ALPHA1", "G_ALPHA2", "C_ALPHA3"],
    "B2M": ["C_B2M"],
    "MHC_II_ALPHA": ["G_ALPHA", "C_IIA"],
    "MHC_II_BETA": ["G_BETA", "C_IIB"],
}
DOMAIN_TO_CHAIN = {d: (ct, i) for ct, ds in CHAIN_TYPES.items() for i, d in enumerate(ds)}
MIN_SCORE = 25.0                   # bits; a domain hit below this is ignored
C_ONLY_MIN = 80.0                  # bits; a chain with C-like hits only (no G domain) needs
                                   # one this strong: TCR / antibody constant domains score
                                   # 20-45 bits, MHC alpha3 / alpha2 / beta2 / b2m >= 120

# Model sets.  MODEL_DIRS (list) = an ensemble: every domain is searched with each set and
# the best-scoring hit is used (with its own set's labels).  Default: the broad set (all
# jawed-vertebrate lineages) + the mammal-focused set; the loop conventions of the
# first set are used for both.
MODEL_DIRS = None
DEFAULT_DIRS = [HERE / "models" / "mhc_broad", HERE / "models" / "mhc_mammal"]

_cache = {}


def _sets():
    """[(hmms, meta)] for every model set in use."""
    if "sets" not in _cache:
        from pyhmmer import plan7
        out = []
        for d in (MODEL_DIRS or ([MODEL_DIR] if MODEL_DIR else DEFAULT_DIRS)):
            d = Path(d)
            with plan7.HMMFile(d / "hmms.hmm") as fh:
                hmms = list(fh)
            out.append((hmms, json.loads((d / "hmms.json").read_text())["models"]))
        # loop conventions: one for all sets (the first set that defines them)
        canon = {}
        for _, meta in out:
            for k, v in meta.items():
                if v.get("loop_canon") and k not in canon:
                    canon[k] = v["loop_canon"]
        _cache["sets"], _cache["canon"] = out, canon
    return _cache["sets"]


def _set_name(k):
    dirs = MODEL_DIRS or ([MODEL_DIR] if MODEL_DIR else DEFAULT_DIRS)
    return Path(dirs[k]).name


def _models(k=0):
    return _sets()[k]


# --------------------------------------------------------------------------- #
#  domain search
# --------------------------------------------------------------------------- #
def find_domains(seq, min_score=MIN_SCORE):
    """All non-overlapping domain hits, best first: list of dicts with model, score,
    and per-residue states: {query_index: ('M', state) | ('I', after_state)}."""
    from pyhmmer import easel, plan7
    hits = []
    for k, (hmms, meta) in enumerate(_sets()):
        abc = hmms[0].alphabet
        dseq = easel.TextSequence(name=b"q", sequence=seq).digitize(abc)
        pipe = plan7.Pipeline(abc, bit_cutoffs=None, E=1e3, domE=1e3)
        for hmm in hmms:
            th = pipe.search_hmm(hmm, easel.DigitalSequenceBlock(abc, [dseq]))
            for h in th:
                for d in h.domains:
                    if d.score < min_score:
                        continue
                    name = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
                    hit = _parse_domain(name, d)
                    hit["set"] = k
                    hits.append(hit)
            pipe.clear()
    hits.sort(key=lambda h: -h["score"])
    chosen = []
    for h in hits:                       # one hit per model region; models may overlap
        span = set(range(h["start"], h["end"] + 1))
        if any(c["model"] == h["model"] and
               len(span & set(range(c["start"], c["end"] + 1))) > 0.1 * len(span)
               for c in chosen):
            continue
        chosen.append(h)
    chosen.sort(key=lambda h: h["start"])
    return chosen


def _parse_domain(model, dom):
    ali = dom.alignment
    hs, ts = ali.hmm_sequence, ali.target_sequence
    state = ali.hmm_from - 1
    qi = ali.target_from - 1
    states = {}
    for hc, tc in zip(hs, ts):
        if hc == ".":                      # insert state: residue after match `state`
            states[qi] = ("I", state)
            qi += 1
        elif tc == "-":                    # delete state
            state += 1
        else:
            state += 1
            states[qi] = ("M", state)
            qi += 1
    return {"model": model, "score": float(dom.score), "start": ali.target_from - 1,
            "end": ali.target_to - 1, "hmm_from": ali.hmm_from, "hmm_to": ali.hmm_to,
            "states": states}


# --------------------------------------------------------------------------- #
#  labels for one domain
# --------------------------------------------------------------------------- #
def _split_label(l):
    i = len(l)
    while i and l[i - 1].isalpha():
        i -= 1
    return int(l[:i]), l[i:]


def _canon_loops(out, windows):
    """Relabel residues between two framing core positions with the canonical IMGT loop
    pattern for that loop length (learned from the IMGT structures; hmms.json)."""
    if not windows:
        return
    by_lab = {lab: qi for qi, (lab, kind) in out.items()}
    for w, by_k in windows.items():
        a, b = (_split_label(x) for x in w.split("|"))
        if a not in by_lab or b not in by_lab:
            continue
        qa, qb = by_lab[a], by_lab[b]
        mid = list(range(qa + 1, qb))
        pat = by_k.get(str(len(mid)))
        if pat is None or not all(q in out for q in mid):
            continue
        for q, l in zip(mid, pat):
            out[q] = (_split_label(l), out[q][1])


def label_domain(hit, n_res):
    """{query_index: ((pos, icode), kind)} for one domain hit; kind 'M' match, 'I' insert,
    'X' terminal extension (residues next to a clipped local alignment, mapped onto the
    unaligned first/last match states of the model)."""
    _, meta = _models(hit.get("set", 0))
    m = meta[hit["model"]]
    match = m["match_labels"]
    M = len(match)
    ins_slots = {int(k): v for k, v in m["insert_labels"].items()}
    out = {}
    runs = {}
    for qi in sorted(hit["states"]):
        kind, st = hit["states"][qi]
        if kind == "M":
            out[qi] = (_split_label(match[st - 1]), "M")
        else:
            runs.setdefault(st, []).append(qi)
    used = set(v[0] for v in out.values())
    for st, qis in runs.items():
        slots = [_split_label(l) for l in ins_slots.get(st, [])]
        slots = [x for x in slots if x not in used]
        prev = out[qis[0] - 1][0] if qis[0] - 1 in out else (
            _split_label(match[st - 1]) if st > 0 else (1, ""))
        for qi, lab in zip(qis, _fill(slots, len(qis), prev)):
            out[qi] = (lab, "I")
            used.add(lab)
    _sets()
    _canon_loops(out, _cache["canon"].get(hit["model"], {}))
    # terminal extension onto unaligned model states
    for k in range(1, hit["hmm_from"]):
        qi, st = hit["start"] - k, hit["hmm_from"] - k
        if qi < 0:
            break
        out[qi] = (_split_label(match[st - 1]), "X")
    for k in range(1, M - hit["hmm_to"] + 1):
        qi, st = hit["end"] + k, hit["hmm_to"] + k
        if qi >= n_res:
            break
        out[qi] = (_split_label(match[st - 1]), "X")
    return out


EXT_MAX = 40          # residues examined beyond a clipped hit end
EXT_GAP = -2.0        # bits, gap in the end-extension DP
EXT_LINK = -0.3      # bits per linker residue between two domains
EXT_MIN_CLIP = 4      # only re-align ends clipped by at least this many model states


def _mk(hit):
    """Model key 'NAME@set' for the emission / transition cache."""
    return f"{hit['model']}@{hit.get('set', 0)}"


def _emission_scores(model):
    """(match log2-odds [states x 20], symbols, log2 transitions [0..M x 7]) of a model.
    Transition columns: MM MI MD IM II DM DD."""
    key = ("emis", model)
    if key not in _cache:
        import numpy as np
        from pyhmmer import plan7
        name, _, k = model.partition("@")
        hmms, _ = _models(int(k or 0))
        model = name
        hmm = next(h for h in hmms if (h.name.decode() if isinstance(h.name, bytes) else h.name) == model)
        E = np.asarray(hmm.match_emissions)[1:, :20]          # rows = states 1..M
        bg = np.asarray(plan7.Background(hmm.alphabet).residue_frequencies)[:20]
        T = np.log2(np.maximum(np.asarray(hmm.transition_probabilities), 1e-6))
        _cache[key] = (np.log2(np.maximum(E, 1e-6) / bg), hmm.alphabet.symbols[:20], T)
    return _cache[key]


def _viterbi(seq, res, states, free_start, free_end, link_after=None, anchor=None):
    """Align residues `res` (chain order) to model states `states` = [(model, k)] (model
    order) with match / insert / delete states and the models' own transitions.
    free_start / free_end: residues and states before the first / after the last aligned
    pair cost nothing (local end); otherwise both lists are anchored there.
    link_after = j: residues between state j and j+1 (a domain junction) are linker 'L'
    at EXT_LINK bits each.  Returns [(residue, 'M'|'I'|'L', (model, k))]."""
    import numpy as np
    nr, ns = len(res), len(states)
    if not nr or not ns:
        return []
    E, T, sym = [], [], None
    for m, k in states:
        Sm, sym, Tm = _emission_scores(m)
        E.append(Sm[k - 1]); T.append(Tm[k])
    idx = {x: i for i, x in enumerate(sym)}
    E, T = np.array(E), np.array(T)
    NEG = -1e9
    # previous state's transition into state j (model boundary: free-ish)
    def tr(j, col):                        # transition out of state j-1 (1-based j>1)
        if j <= 1:
            return 0.0
        if states[j - 2][0] != states[j - 1][0]:
            return 0.0                     # crossing into the next domain
        return T[j - 2, col]
    Mx = np.full((nr + 1, ns + 1), NEG); Ix = np.full((nr + 1, ns + 1), NEG)
    Dx = np.full((nr + 1, ns + 1), NEG)
    bt = {}
    Ta = None
    if anchor is not None:                 # transitions out of the anchor state (hit end)
        Ta = _emission_scores(anchor[0])[2][anchor[1]]
    Mx[0, 0] = 0.0
    if free_start:
        Mx[0, :] = 0.0; Mx[:, 0] = 0.0
    else:
        Ix[:, 0] = [EXT_LINK * i if link_after == 0 else EXT_GAP * i for i in range(nr + 1)]
        if Ta is not None:                 # MI then II from the anchor
            Ix[1:, 0] = Ta[1] + Ta[4] * np.arange(nr)
        Ix[0, 0] = NEG
    for i in range(0, nr + 1):
        a = idx.get(seq[res[i - 1]]) if i else None
        for j in range(0, ns + 1):
            if i == 0 and j == 0:
                continue
            if i and j:
                e = E[j - 1, a] if a is not None else 0.0
                c = [(Mx[i - 1, j - 1] + tr(j, 0), "M"), (Ix[i - 1, j - 1] + tr(j, 3), "I"),
                     (Dx[i - 1, j - 1] + tr(j, 5), "D")]
                if j == 1 and not free_start:
                    c.append((Ix[i - 1, 0] + (Ta[3] if Ta is not None else 0.0), "I"))
                v, f = max(c)
                if free_start and i - 1 >= 0 and j == 1:
                    pass
                if v + e > Mx[i, j]:
                    Mx[i, j] = v + e; bt["M", i, j] = f
            if i and j:
                link = link_after == j
                c = [(Mx[i - 1, j] + (EXT_LINK if link else T[j - 1, 1]), "M"),
                     (Ix[i - 1, j] + (EXT_LINK if link else T[j - 1, 4]), "I")]
                v, f = max(c)
                Ix[i, j] = v; bt["I", i, j] = f
            if j:
                c = [(Mx[i, j - 1] + tr(j, 2), "M"), (Dx[i, j - 1] + tr(j, 6), "D")]
                if j == 1:
                    # entering the state list by a deletion: from the anchor or after
                    # residues inserted next to the anchor (e.g. carp alpha3 K2111, then
                    # 2113-2115 deleted)
                    if Ta is not None:
                        c = [(Mx[i, 0] + Ta[2], "M"), (Ix[i, 0] + Ta[3] + Ta[2], "I")]
                    else:
                        c = [(Mx[i, 0] + (0.0 if free_start else EXT_GAP), "M"),
                             (Ix[i, 0] + EXT_GAP, "I")]
                v, f = max(c)
                Dx[i, j] = v; bt["D", i, j] = f
    if free_end:
        cands = [(Mx[i, j], "M", i, j) for i in range(nr + 1) for j in range(ns + 1)]
        v, cur, i, j = max(cands)
    else:
        v, cur = max((Mx[nr, ns], "M"), (Ix[nr, ns], "I"), (Dx[nr, ns], "D"))
        i, j = nr, ns
    out = []
    while i > 0 or j > 0:
        if cur == "M":
            if j == 0:                       # residues before any state: free start or linker
                break
            if i == 0:
                break
            out.append((res[i - 1], "M", states[j - 1]))
            prev = bt.get(("M", i, j))
            if prev is None:
                break
            i -= 1; j -= 1; cur = prev
        elif cur == "I":
            if j == 0:      # next to the anchor: linker at a junction, else an insertion
                out.append((res[i - 1], "L", None) if link_after is not None
                           else (res[i - 1], "I", states[0])); i -= 1
                if i == 0:
                    break
                continue
            kind = "L" if link_after == j else "I"
            out.append((res[i - 1], kind, states[j - 1]))
            cur = bt["I", i, j]; i -= 1
        else:
            prev = bt.get(("D", i, j), "M")
            j -= 1; cur = prev
            if j == 0 and cur == "M":
                break
    return out[::-1]


def _apply(hit, pairs):
    """Add aligned residues [(residue, 'M'|'I', state)] to a hit."""
    st = hit["states"]
    for qi, kind, s_ in pairs:
        if kind == "M":
            st[qi] = ("M", s_)
    for qi, kind, s_ in pairs:
        if kind == "I":
            prev = [st[q][1] for q in st if q < qi and st[q][0] == "M"]
            st[qi] = ("I", max(prev) if prev else 0)
    hit["start"], hit["end"] = min(st), max(st)
    hit["hmm_from"] = min(v[1] for v in st.values() if v[0] == "M")
    hit["hmm_to"] = max(v[1] for v in st.values() if v[0] == "M")


def _extend_ends(hit, seq, lo, hi):
    """Re-align residues lo..start-1 to states 1..hmm_from-1 (anchored at the hit, free
    N-terminal end) or end+1..hi to hmm_to+1..M (anchored at the hit, free C-terminal
    end): local alignments are often clipped at domain ends."""
    m = _mk(hit)
    M = _emission_scores(m)[0].shape[0]
    pairs = []
    if hit["hmm_from"] - 1 >= EXT_MIN_CLIP and lo < hit["start"]:
        pairs += _viterbi(seq, list(range(lo, hit["start"])),
                          [(m, k) for k in range(1, hit["hmm_from"])], True, False)
    if M - hit["hmm_to"] >= EXT_MIN_CLIP and hi > hit["end"]:
        pairs += _viterbi(seq, list(range(hit["end"] + 1, hi + 1)),
                          [(m, k) for k in range(hit["hmm_to"] + 1, M + 1)], False, True,
                          anchor=(m, hit["hmm_to"]))
    pairs = [(q, kd, s[1]) for q, kd, s in pairs if kd in "MI"]
    if pairs:
        _apply(hit, pairs)
    return hit


def _extend_junction(a, b, seq):
    """Share the residues between hit a and hit b (a before b) over a's unaligned tail
    states and b's unaligned head states: one alignment anchored at both hits, linker
    residues allowed only at the junction."""
    res = list(range(a["end"] + 1, b["start"]))
    if not res:
        return
    Ma = _emission_scores(_mk(a))[0].shape[0]
    tail = list(range(a["hmm_to"] + 1, Ma + 1)) if Ma - a["hmm_to"] >= EXT_MIN_CLIP else []
    head = list(range(1, b["hmm_from"])) if b["hmm_from"] - 1 >= EXT_MIN_CLIP else []
    if not tail and not head:
        return
    states = [(_mk(a), k) for k in tail] + [(_mk(b), k) for k in head]
    out = _viterbi(seq, res, states, False, False, link_after=len(tail))
    pa = [(q, kd, s[1]) for q, kd, s in out if kd in "MI" and s[0] == _mk(a) and
          (s[1] in tail) and states.index(s) < len(tail)]
    pb = [(q, kd, s[1]) for q, kd, s in out if kd in "MI" and states.index(s) >= len(tail)]
    if pa:
        _apply(a, pa)
    if pb:
        _apply(b, pb)


def _priority(model, label, kind):
    """Which domain owns a residue claimed by two domain hits."""
    pos, ic = label
    r = {"M": 2.0, "I": 1.0, "X": 0.0}[kind]
    if model.startswith("G_") and pos > 92:          # G grid tail (93-97): loses to next domain
        r -= 3.0
    if pos == 1 and ic:                              # N-terminal insertion positions: weak
        r -= 0.5
    return r


def _fill(slots, n, prev):
    """Choose n labels from the defined insertion slots (in column order)."""
    if n <= len(slots):
        positions = sorted(set(p for p, _ in slots))
        if len(positions) == 2 and positions[1] == positions[0] + 1 and all(
                ic for _, ic in slots):
            # symmetric loop (C-LIKE 84A.. / ..85A): outer letters first on both sides
            left = [s for s in slots if s[0] == positions[0]]
            right = [s for s in slots if s[0] == positions[1]]     # listed ..85B, 85A
            n_r = n // 2
            n_l = n - n_r
            if n_l <= len(left) and n_r <= len(right):
                return left[:n_l] + (right[len(right) - n_r:] if n_r else [])
        return slots[:n]
    # more residues than defined slots: use all slots, then extra letters after the last one
    labels = list(slots)
    base = labels[-1] if labels else prev
    p, ic = base
    k = LETTERS.index(ic) + 1 if ic else 0
    while len(labels) < n:
        labels.append((p, LETTERS[k]))
        k += 1
    return labels


# --------------------------------------------------------------------------- #
#  chain numbering
# --------------------------------------------------------------------------- #
def classify(hits):
    """Chain type from the domain hits (by summed score); None if no MHC chain (no hit, or
    only weak C-like hits such as a TCR or antibody constant domain)."""
    if not any(h["model"].startswith("G_") for h in hits) and \
            max([h["score"] for h in hits] + [0]) < C_ONLY_MIN:
        return None
    score = {}
    for h in hits:
        ct, _ = DOMAIN_TO_CHAIN[h["model"]]
        score[ct] = score.get(ct, 0) + h["score"]
    return max(score, key=score.get) if score else None


def number_sequence(seq, chain_type=None, gaps=None):
    """IMGT-number an MHC chain sequence.

    gaps: {i: n} = n residues are missing (chain break) before residue i. They are
    represented by 'X' placeholders while aligning, so the HMM puts the gap where the
    coordinates are missing, and removed again afterwards.

    Returns (resids, info): resids = [(resseq:int, icode:str)] for every residue (unique),
    info = dict(chain_type, domains=[(model, start, end, score)], extra_insertions)."""
    if gaps:
        return _number_with_gaps(seq, chain_type, gaps)
    hits = find_domains(seq)
    ct = chain_type or classify(hits)
    if ct is None:
        raise ValueError("no MHC domain found")
    hits = [h for h in hits if DOMAIN_TO_CHAIN[h["model"]][0] == ct]
    # one hit per domain model (tandem copies, e.g. duplicated class I genes in one ORF):
    # best-scoring hits that are in domain order; the other copy is numbered as tail
    best = {}
    for h in sorted(hits, key=lambda h: -h["score"]):
        best.setdefault(h["model"], h)
    kept = []
    for h in sorted(best.values(), key=lambda h: -h["score"]):
        d = DOMAIN_TO_CHAIN[h["model"]][1]
        if all((DOMAIN_TO_CHAIN[k["model"]][1] < d) == (k["start"] < h["start"]) for k in kept):
            kept.append(h)
    hits = sorted(kept, key=lambda h: h["start"])
    n = len(seq)
    # local alignments are often clipped at the domain ends; extend them by aligning the
    # neighbouring residues to the unaligned model states (emission scores, anchored DP)
    # residues between two hits are shared out by ONE joint alignment (tail states of the
    # first domain, then head states of the second), so neither domain eats the other's ends
    if hits:
        _extend_ends(hits[0], seq, max(0, hits[0]["start"] - EXT_MAX), hits[0]["start"] - 1)
        _extend_ends(hits[-1], seq, hits[-1]["end"] + 1, min(n - 1, hits[-1]["end"] + EXT_MAX))
    for a, b in zip(hits, hits[1:]):
        _extend_junction(a, b, seq)
    lab = [None] * n
    best = [None] * n                                # (priority, -domain index)
    for h in sorted(hits, key=lambda h: DOMAIN_TO_CHAIN[h["model"]][1]):
        d_idx = DOMAIN_TO_CHAIN[h["model"]][1]
        for qi, ((p, ic), kind) in label_domain(h, n).items():
            pr = (_priority(h["model"], (p, ic), kind), -d_idx)
            if best[qi] is None or pr > best[qi]:
                best[qi] = pr
                lab[qi] = (1000 * d_idx + p, ic)
    # labels must increase along the chain domain-wise; drop claims that go backwards
    for i in range(1, n):
        if lab[i] is not None:
            prev = next((lab[j] for j in range(i - 1, -1, -1) if lab[j] is not None), None)
            if prev is not None and lab[i][0] // 1000 < prev[0] // 1000:
                lab[i] = None
    seen = set()
    for i in range(n):                               # duplicates (should be rare): keep first
        if lab[i] is not None:
            if lab[i] in seen:
                lab[i] = None
            else:
                seen.add(lab[i])
    numbered = [i for i, l in enumerate(lab) if l is not None]
    if not numbered:
        raise ValueError("no residue numbered")
    first, last = numbered[0], numbered[-1]
    used = set(l for l in lab if l)
    extra = 0
    i = first + 1
    while i < last:                              # linkers / unaligned residues inside
        if lab[i] is not None:
            i += 1
            continue
        j = i
        while lab[j] is None:
            j += 1
        p, ic = lab[i - 1]
        k = LETTERS.index(ic) + 1 if ic else 0
        free = [L for L in LETTERS[k:] if (p, L) not in used]
        if j - i <= min(10, len(free)):
            # short run: insertion codes on the preceding position
            for q, L in zip(range(i, j), free):
                lab[q] = (p, L)
        else:
            # long unaligned segment (not part of any MHC domain): linker block 900+ of the
            # preceding domain, e.g. 1901, 1902, ...
            base = 1000 * ((p - 1) // 1000 if p > 0 else 0) + 900
            n0 = 1 + max([x[0] - base for x in used if base < x[0] < base + 100] or [0])
            for q in range(i, j):
                lab[q] = (base + n0 + (q - i), "")
        for q in range(i, j):
            used.add(lab[q])
        extra += j - i
        i = j
    for k, i in enumerate(range(first - 1, -1, -1)):
        lab[i] = (-k, "")
    tail = 1000 * len(CHAIN_TYPES[ct])
    for k, i in enumerate(range(last + 1, n), start=1):
        lab[i] = (tail + k, "")
    if len(set(lab)) != n:
        # safety net (unusual chains, e.g. internal repeats): later duplicates go to the
        # 950+ block of their domain, reported in info["renumbered_duplicates"]
        seen, dup = set(), 0
        for i, l in enumerate(lab):
            if l in seen:
                base = 1000 * ((l[0] - 1) // 1000 if l[0] > 0 else 0) + 950
                k = 1
                while (base + k, "") in seen or (base + k, "") in lab:
                    k += 1
                lab[i] = (base + k, "")
                dup += 1
            seen.add(lab[i])
        extra_info = {"renumbered_duplicates": dup}
    else:
        extra_info = {}
    info = {"chain_type": ct, "extra_insertions": extra, **extra_info,
            "domains": [(h["model"], h["start"], h["end"], round(h["score"], 1)) for h in hits],
            "model_sets": [_set_name(h.get("set", 0)) for h in hits]}
    return lab, info


PAD_CHOICES = (1, 2, 3, 4, 6, 8, 10, 12, 15, 20)


def _padded(seq, pads):
    padded, keep = [], []
    for i, c in enumerate(seq):
        if pads.get(i):
            padded += ["X"] * pads[i]
            keep += [False] * pads[i]
        padded.append(c)
        keep.append(True)
    return "".join(padded), keep


def _total_score(seq, chain_type):
    hits = find_domains(seq)
    ct = chain_type or classify(hits)
    return sum(h["score"] for h in hits if DOMAIN_TO_CHAIN[h["model"]][0] == ct), ct


def _number_with_gaps(seq, chain_type, gaps):
    """Missing residues at chain breaks are represented by 'X' placeholders. How many is
    chosen per break (coordinate ascent over PAD_CHOICES, starting from the estimate from
    the CA-CA distance) by the total HMM score: the models know how long each loop is."""
    pads = {i: max(1, n) for i, n in gaps.items()}
    best, ct = _total_score(_padded(seq, pads)[0], chain_type)
    if ct is None:
        raise ValueError("no MHC domain found")
    for _ in range(2):
        changed = False
        for i in sorted(pads):
            for n in PAD_CHOICES:
                if n < gaps[i] or n == pads[i]:      # never fewer than the distance allows
                    continue
                trial = dict(pads)
                trial[i] = n
                sc, _ = _total_score(_padded(seq, trial)[0], ct)
                if sc > best + 0.5:
                    best, pads, changed = sc, trial, True
        if not changed:
            break
    padded, keep = _padded(seq, pads)
    lab, info = number_sequence(padded, ct)
    info["gap_fill"] = {int(k): v for k, v in pads.items()}
    return [l for l, k in zip(lab, keep) if k], info


def chain_gaps(residues, ca_step=3.8, max_fill=30):
    """{i: n_missing} for a list of Bio.PDB residues: chain breaks (C-N > 2 A) with the
    number of missing residues estimated from the CA-CA distance (at least 1)."""
    import math
    import numpy as np
    gaps = {}
    for i in range(1, len(residues)):
        a, b = residues[i - 1], residues[i]
        if "C" in a and "N" in b and np.linalg.norm(a["C"].coord - b["N"].coord) <= 2.0:
            continue
        if "CA" in a and "CA" in b:
            d = float(np.linalg.norm(a["CA"].coord - b["CA"].coord))
            if d < 4.2:
                continue
            gaps[i] = min(max_fill, max(1, math.ceil(d / ca_step) - 1))
    return gaps


def number_residues(residues, chain_type=None):
    """IMGT numbering for a list of Bio.PDB residues (gap-aware)."""
    from Bio.Data.PDBData import protein_letters_3to1
    seq = "".join(protein_letters_3to1.get(r.get_resname(), "X") for r in residues)
    return number_sequence(seq, chain_type, chain_gaps(residues))
