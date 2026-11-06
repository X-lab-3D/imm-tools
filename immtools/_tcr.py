# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 <Li Xue>
"""TCR (and antibody) variable-domain numbering with ANARCI (package anarci-pyhmmer: the
X-lab-3D fork of David Prihoda's pyhmmer port, IMGT HMMs included).
"""
import os
import sys

CHAIN_TYPE_LABELS = {
    "A": "TCR alpha (Valpha)",
    "B": "TCR beta (Vbeta)",
    "G": "TCR gamma (Vgamma)",
    "D": "TCR delta (Vdelta)",
    "H": "antibody heavy (VH)",
    "K": "antibody kappa (VL-kappa)",
    "L": "antibody lambda (VL-lambda)",
}


def available():
    try:
        from anarci import anarci  # noqa: F401
        return True
    except ImportError:
        return False


def status():
    """(ok, message) about the ANARCI installation."""
    import os
    try:
        import anarci as pkg
        from anarci import anarci  # noqa: F401
    except ImportError as e:
        return False, f"ANARCI not importable ({e}); install it with: pip install anarci-pyhmmer"
    hmm = os.path.join(os.path.dirname(pkg.__file__), "dat", "HMMs", "ALL.hmm")
    if not os.path.exists(hmm):
        return False, (f"ANARCI found at {os.path.dirname(pkg.__file__)} but without its HMMs "
                       "(another anarci package?); reinstall with: pip install --force-reinstall anarci-pyhmmer")
    r = number_chain("DAGVIQSPRHEVTEMGQEVTLRCKPISGHNSLFWYRQTMMRGLELLIYFNNNVPIDDSGMPEDRFSAKMPNASFSTLKIQPSEPRDSAVYFCASSYVGNTGELFFGEGSRLTVL")
    if not r or r[1][0]["chain_type"] != "B":
        return False, "ANARCI is installed but did not number a test TCR beta chain"
    return True, f"ANARCI OK ({os.path.dirname(pkg.__file__)})"


def require():
    if not available():
        sys.stderr.write("ERROR!! TCR numbering needs ANARCI:  pip install anarci-pyhmmer\n")
        sys.exit(1)


def number_chain(seq, scheme="imgt", bit_score_threshold=80):
    """ANARCI on one chain sequence.

    Returns None (no variable domain) or (labels, domains): labels = [(number, icode) or None]
    per residue of `seq` (None = outside the variable domain); domains = [dict(chain_type,
    label, species, start, end, bitscore)]."""
    from anarci import anarci
    numbering, details, _ = anarci([("q", seq)], scheme=scheme, output=False,
                                   bit_score_threshold=bit_score_threshold,
                                   assign_germline=True)
    if numbering[0] is None:
        return None
    labels = [None] * len(seq)
    doms = []
    pairs = sorted(zip(numbering[0], details[0]), key=lambda x: x[1]["query_start"])
    for k, (dom_num, det) in enumerate(pairs):
        # several V domains in one chain (single-chain constructs): 2nd one +1000, ...
        i = det["query_start"]
        for (num, icode), aa in dom_num[0]:
            if aa == "-":
                continue
            labels[i] = (1000 * k + num, icode.strip() if icode else "")
            i += 1
        g = det.get("germlines") or {}
        doms.append(dict(chain_type=det["chain_type"],
                         v_gene=(g.get("v_gene") or [[None, None]])[0][1],
                         v_identity=(g.get("v_gene") or [None, None])[1],
                         j_gene=(g.get("j_gene") or [[None, None]])[0][1],
                         j_identity=(g.get("j_gene") or [None, None])[1],
                         label=CHAIN_TYPE_LABELS.get(det["chain_type"], det["chain_type"]),
                         species=det["species"], start=det["query_start"],
                         end=det["query_end"] - 1, bitscore=det.get("bitscore")))
    return labels, doms


def continue_after(labels):
    """Residues outside the V domain (constant domain, tags): sequential after the highest
    IMGT number.  Residues before the V domain get
    0, -1, ... counting backwards."""
    out = list(labels)
    nums = [l[0] for l in out if l]
    if not nums:
        return out
    first = next(i for i, l in enumerate(out) if l)
    for k, i in enumerate(range(first - 1, -1, -1)):
        out[i] = (-k, "")
    nxt = max(nums) + 1
    for i in range(first, len(out)):
        if out[i] is None:
            out[i] = (nxt, "")
            nxt += 1
    return out


# --------------------------------------------------------------------------- #
#  constant domain (IMGT C-LIKE numbering, profile HMMs built from IMGT structures)
# --------------------------------------------------------------------------- #
_C = {}
C_MIN_SCORE = 40.0


def _c_models():
    if "hmms" not in _C:
        import json
        from pathlib import Path
        from pyhmmer import plan7
        d = Path(__file__).resolve().parent / "models" / "tcr_constant"
        with plan7.HMMFile(d / "mhc_hmms.hmm") as fh:
            _C["hmms"] = list(fh)
        _C["meta"] = json.loads((d / "mhc_hmms.json").read_text())["models"]
    return _C["hmms"], _C["meta"]


def number_constant(seq, start):
    """IMGT C-LIKE labels for the constant domain found in seq[start:].
    Returns ({index: (position, icode)}, model, score) or (None, None, None)."""
    from pyhmmer import easel, plan7
    from . import _mhc
    hmms, meta = _c_models()
    abc = hmms[0].alphabet
    sub = seq[start:]
    if len(sub) < 40:
        return None, None, None
    dseq = easel.TextSequence(name=b"q", sequence=sub).digitize(abc)
    best = None
    pipe = plan7.Pipeline(abc, bit_cutoffs=None, E=1e3, domE=1e3)
    for hmm in hmms:
        name = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
        for h in pipe.search_hmm(hmm, easel.DigitalSequenceBlock(abc, [dseq])):
            for d in h.domains:
                if d.score >= C_MIN_SCORE and (best is None or d.score > best[1].score):
                    best = (name, d)
        pipe.clear()
    if best is None:
        return None, None, None
    name, d = best
    hit = _mhc._parse_domain(name, d)
    m = meta[name]
    match = m["match_labels"]
    slots = {int(k): v for k, v in m["insert_labels"].items()}
    out, runs = {}, {}
    for qi in sorted(hit["states"]):
        kind, s = hit["states"][qi]
        if kind == "M":
            out[qi] = _mhc._split_label(match[s - 1])
        else:
            runs.setdefault(s, []).append(qi)
    used = set(out.values())
    for s, qis in runs.items():
        free = [_mhc._split_label(x) for x in slots.get(s, [])]
        free = [x for x in free if x not in used]
        prev = out.get(qis[0] - 1) or (_mhc._split_label(match[s - 1]) if s else (1, ""))
        for qi, lab in zip(qis, _mhc._fill(free, len(qis), prev)):
            out[qi] = lab
            used.add(lab)
    # clipped ends: unaligned first / last match states, contiguous
    for k in range(1, hit["hmm_from"]):
        qi = hit["start"] - k
        if qi < 0:
            break
        out[qi] = _mhc._split_label(match[hit["hmm_from"] - k - 1])
    for k in range(1, len(match) - hit["hmm_to"] + 1):
        qi = hit["end"] + k
        if qi >= len(sub):
            break
        out[qi] = _mhc._split_label(match[hit["hmm_to"] + k - 1])
    return {start + q: l for q, l in out.items()}, name, round(float(d.score), 1)


def full_labels(seq, constant=False):
    """Labels for every residue of a TCR chain, or None if no variable domain.

    V domain: IMGT (ANARCI).  constant=False: the rest continues after the highest IMGT
    number (residues before the V domain: 0, -1, ...).  constant=True: the constant domain
    gets IMGT C-LIKE positions + 1000 (1001-1128, as IMGT/3Dstructure-DB), residues between
    V and C domain continue the V domain, residues after the C domain 2001, 2002, ...
    Returns (labels, info) with info = dict(domains=[...], constant=(model, score) | None)."""
    r = number_chain(seq)
    if r is None:
        return None
    vlab, doms = r
    info = dict(domains=doms, constant=None)
    if not constant:
        return continue_after(vlab), info
    vend = max(i for i, l in enumerate(vlab) if l)
    cblock = 1000 * (1 + max(l[0] for l in vlab if l) // 1000)    # next free thousand
    clab, model, score = number_constant(seq, vend + 1)
    if not clab:
        return continue_after(vlab), info
    info["constant"] = (model, score)
    out = list(vlab)
    for i, (p, ic) in clab.items():
        if out[i] is None:
            out[i] = (cblock + p, ic)
    cend = max(clab)
    first = next(i for i, l in enumerate(out) if l)
    for k, i in enumerate(range(first - 1, -1, -1)):
        out[i] = (-k, "")
    # V-C linker residues not labelled: continue V with letters; after C: 2001...
    nxt = cblock + 1001
    for i in range(first, len(out)):
        if out[i] is None:
            if i > cend:
                out[i] = (nxt, ""); nxt += 1
            else:
                prev = out[i - 1]
                letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                ic = letters[letters.index(prev[1]) + 1] if prev[1] else "A"
                out[i] = (prev[0], ic)
    # uniqueness (should hold); fall back to V-only numbering if not
    if len(set(out)) != len(out):
        return continue_after(vlab), dict(domains=doms, constant=None)
    return out, info
