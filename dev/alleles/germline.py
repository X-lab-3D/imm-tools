"""Shared helpers for the allele-frequency source builders (dev/alleles/build_*.py).

A builder reads one external resource and calls write_source() with metadata and rows; the
table format is described in immtools/_alleles.py. The helpers here turn a germline protein
into the 128-column IMGT string by aligning it to ANARCI's IMGT-gapped human germline of the
same gene, so every source is numbered the same way.
"""
import datetime
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from immtools._alleles import COLUMNS, DIR, POSITIONS, gene_key   # noqa: E402

META_KEYS = ("name", "title", "species", "populations", "citation", "url", "licence")
_TEMPLATES = {}


def templates(species="human"):
    """segment -> gene -> ANARCI's IMGT-gapped germline (first allele of the gene)."""
    if species not in _TEMPLATES:
        from anarci.germlines import all_germlines
        t = {"V": {}, "J": {}}
        for seg in "VJ":
            for chain in "ABGD":
                for name, g in all_germlines[seg].get(chain, {}).get(species, {}).items():
                    t[seg].setdefault(gene_key(name), g)
        _TEMPLATES[species] = t
    return _TEMPLATES[species]


def _aligner():
    from Bio.Align import PairwiseAligner, substitution_matrices
    al = PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-1,
                         substitution_matrix=substitution_matrices.load("BLOSUM62"))
    al.end_gap_score = 0              # germline sequences may be partial at either end
    return al


def number(protein, gene, segment, species="human", min_identity=0.7):
    """128-column IMGT string of the germline-encoded part of `protein`, or None if ANARCI has
    no germline of `gene` or the protein aligns poorly to it."""
    g = templates(species)[segment].get(gene_key(gene))
    if not g:
        return None
    cols = [i for i, c in enumerate(g) if c != "-"]
    t = "".join(g[i] for i in cols)
    aln = _aligner().align(protein, t)[0]
    keep = set(POSITIONS[segment])
    out, same, n = ["-"] * 128, 0, 0
    for (q0, q1), (t0, _) in zip(*aln.aligned):
        for k in range(q1 - q0):
            pos = cols[t0 + k] + 1
            if pos in keep:
                out[pos - 1] = protein[q0 + k]
                n += 1
                same += protein[q0 + k] == t[t0 + k]
    if n < (40 if segment == "V" else 6) or same / n < min_identity:
        return None
    return "".join(out)


def translations(protein, nucleotides=""):
    """Candidate proteins: the given one if it has no stop codon, else the reading frames of
    `nucleotides` without one (a leading stop from a partial first codon is dropped).
    Empty: a stop in every frame, i.e. not a functional allele."""
    protein = protein.lstrip("*")
    if "*" not in protein:
        return [protein]
    from Bio.Seq import Seq
    out = []
    for f in range(3):
        s = nucleotides[f:]
        p = str(Seq(s[:len(s) // 3 * 3]).translate()).lstrip("*")
        if p and "*" not in p:
            out.append(p)
    return out


def write_source(meta, rows, script):
    """Write DIR/<meta['name']>.tsv. rows: dicts with the keys of COLUMNS."""
    missing = [k for k in META_KEYS if not meta.get(k)]
    if missing:
        raise ValueError(f"metadata missing: {missing}")
    path = DIR / f"{meta['name']}.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted({tuple(str(r[k]) for k in COLUMNS) for r in rows})
    with open(path, "w") as fh:
        for k in META_KEYS:
            fh.write(f"# {k}: {meta[k]}\n")
        fh.write(f"# built: {datetime.date.today()} by dev/alleles/{Path(script).name}\n")
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return path, len(rows)
