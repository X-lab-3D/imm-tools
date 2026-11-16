# SPDX-License-Identifier: Apache-2.0
"""Closest MHC allele / gene: the chain's peptide-binding region (class I alpha1+alpha2,
class II alpha1 or beta1; whole domain for beta2m) against reference sequences from
IPD-IMGT/HLA, IPD-MHC and UniProt (models/mhc_refs.fasta.gz)."""
import gzip
from pathlib import Path

REFS = Path(__file__).resolve().parent / "models" / "mhc_refs.fasta.gz"
REGION = {"MHC_I_ALPHA": (0, 1), "MHC_II_ALPHA": (0, 0), "MHC_II_BETA": (0, 0), "B2M": (0, 0)}
_cache = {}


K = 3


def _kmers(s):
    return {s[i:i + K] for i in range(len(s) - K + 1)}


def _db(chain_type):
    """(names, sequences, inverted 3-mer index) for one chain type."""
    if chain_type not in _cache:
        import collections
        import numpy as np
        seqs, names, name = [], [], None
        with gzip.open(REFS, "rt") as fh:
            for ln in fh:
                ln = ln.rstrip("\n")
                if ln.startswith(">"):
                    name = ln[1:]
                elif name and name.split("|")[0] == chain_type:
                    names.append(name)
                    seqs.append(ln)
                    name = None
        index = collections.defaultdict(list)
        for i, sq in enumerate(seqs):
            for k in _kmers(sq):
                index[k].append(i)
        index = {k: np.array(v, dtype=np.int32) for k, v in index.items()}
        _cache[chain_type] = (names, seqs, index)
    return _cache[chain_type]


def _identity(q, t):
    """Identical residues / query length, from a global alignment (free end gaps)."""
    from Bio.Align import PairwiseAligner, substitution_matrices
    if "aligner" not in _cache:
        a = PairwiseAligner()
        a.mode = "global"
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score, a.extend_gap_score = -10, -0.5
        a.end_gap_score = 0                       # works in old and new Biopython
        _cache["aligner"] = a
    aln = _cache["aligner"].align(q, t)[0]
    same = sum(q[i] == t[j] for (a0, a1), (b0, b1) in zip(*aln.aligned)
               for i, j in zip(range(a0, a1), range(b0, b1)))
    return same / len(q)


def region(seq, labels, chain_type):
    """Residues of the reference region, from the MHC numbering (domain index = number//1000)."""
    d0, d1 = REGION[chain_type]
    return "".join(a for a, l in zip(seq, labels) if l and d0 <= l[0] // 1000 <= d1
                   and 0 < l[0] % 1000 < 100 + 30 * (chain_type == "B2M"))


def closest(seq, labels, chain_type, shortlist=40):
    """(name, species, identity %, region length) of the closest reference, or None.
    3-mer shortlist, then global alignment identity over the query region."""
    import numpy as np
    q = region(seq, labels, chain_type)
    if len(q) < 40:
        return None
    names, seqs, index = _db(chain_type)
    if not seqs:
        return None
    counts = np.zeros(len(seqs), dtype=np.int32)
    for k in _kmers(q):
        if k in index:
            counts[index[k]] += 1
    cand = np.argsort(-counts)[:shortlist]
    best = max(((_identity(q, seqs[i]), -i) for i in cand), default=None)
    if best is None:
        return None
    _, label, species = names[-best[1]].split("|", 2)
    import re
    label = re.sub(r"\s*\[Cleaved into:[^\]]*\]", "", label)
    return label, species, round(100 * best[0], 1), len(q)
