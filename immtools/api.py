# SPDX-License-Identifier: Apache-2.0
"""Python interface of imm-tools.

Everything the command-line tools do, callable from Python. Structures can be given as a
file path (PDB or mmCIF), an open file, or the file content as a string; functions that
change a structure return the new file content (same format as the input) and write it to
`out` if given.

    import immtools

    immtools.number_mhc("GSHSMRYFFTSVSRPGRGEPRFIAVGYVDD...")   # sequence -> Numbering
    immtools.number_tcr("DAGVIQSPRHEVTEMGQEVTLRC...", constant=True)

    immtools.annotate("3qdg.pdb")        # [{'chain': 'A', 'type': 'MHC_I', ...}, ...]
    immtools.complexes("3qdg.pdb")       # [{'mhc': 'A', 'b2m': 'B', 'pep': 'C', ...}]
    immtools.cdrs("3qdg.pdb")            # [{'chain': 'D', 'cdr': 'CDR3', 'sequence': ...}]
    immtools.check("3qdg_imgt.pdb")      # [{'chain': 'A', 'check': 'anchors', 'status': ...}]

    text = immtools.renumber("3qdg.pdb", out="3qdg_imgt.pdb", constant=True)
    text = immtools.rechain("3qdg.pdb", out="3qdg_std.pdb", layout={"b2m": "L"})
    text = immtools.select_complex("1d9k.pdb", 2)

Stability: the names exported from `immtools` (this module) are the public interface;
modules starting with an underscore are internal and may change.
"""
import io
import os
import sys
from contextlib import redirect_stderr
from dataclasses import dataclass, field
from pathlib import Path

from . import _doc, _io, _mhc, _tcr, _chains, _complex, _number, _registry

__all__ = ["Numbering", "number_mhc", "number_tcr", "annotate", "complexes", "cdrs", "check",
           "renumber", "rechain", "select_complex", "structure"]


@dataclass
class Numbering:
    """IMGT numbering of one chain.

    chain_type: MHC_I, MHC_II_alpha, MHC_II_beta, B2M, TCR_alpha, TCR_beta, ...
    numbers:    (number, insertion code) for every residue of `sequence`
    domains:    [(domain name, first index, last index, bit score)] (0-based indices)
    genes:      TCR only: V / J genes, species (from ANARCI)
    """
    chain_type: str
    sequence: str
    numbers: list
    domains: list = field(default_factory=list)
    genes: dict = field(default_factory=dict)

    @property
    def labels(self):
        """Numbers as strings, e.g. '1001A'."""
        return [f"{n}{ic}" for n, ic in self.numbers]

    def __iter__(self):
        """(residue, label) pairs."""
        return iter(zip(self.sequence, self.labels))


# ------------------------------------------------------------------ sequences

def number_mhc(sequence, chain_type=None):
    """IMGT numbering of an MHC chain (class I heavy chain, class II alpha / beta) or beta2m.

    chain_type: force 'MHC_I_ALPHA', 'MHC_II_ALPHA', 'MHC_II_BETA' or 'B2M'.
    Raises ValueError if no MHC domain is found."""
    seq = sequence.strip().upper()
    labels, info = _mhc.number_sequence(seq, chain_type)
    return Numbering(_registry.MHC_KIND[info["chain_type"]], seq, list(labels),
                     [tuple(d) for d in info["domains"]])


def number_tcr(sequence, constant=False):
    """IMGT numbering of a TCR chain: V domain (ANARCI), and with constant=True the constant
    domain (1001-1128). Raises ValueError if no TCR variable domain is found."""
    _require_anarci()
    seq = sequence.strip().upper()
    r = _tcr.full_labels(seq, constant)
    if r is None:
        raise ValueError("no TCR variable domain found")
    labels, info = r
    d = info["domains"][0]
    kind = _registry.ANARCI_KIND.get(d["chain_type"], "other")
    if _registry.family(kind) != "tcr":
        raise ValueError(f"not a TCR chain ({kind}); antibodies are not numbered yet")
    score = lambda x: None if x.get("bitscore") is None else round(float(x["bitscore"]), 1)
    doms = [("V", x["start"], x["end"], score(x)) for x in info["domains"]]
    if info["constant"]:
        ci = [i for i, l in enumerate(labels) if 1000 < l[0] < 2000]
        if ci:
            doms.append((info["constant"][0], ci[0], ci[-1], info["constant"][1]))
    genes = {k: d.get(k) for k in ("species", "v_gene", "v_identity", "j_gene", "j_identity")}
    return Numbering(kind, seq, list(labels), doms, genes)


# ------------------------------------------------------------------ structures

def annotate(structure, chains=None):
    """imm_annotate as a list of dicts (one per chain)."""
    from . import imm_annotate
    doc, name = _read(structure)
    return list(imm_annotate.rows(doc.model(), _chainset(chains), name))


def complexes(structure):
    """TCR-pMHC complexes: list of dicts role -> chain ID (roles tcra, tcrb, mhc, mhcb, b2m,
    pep), in the order imm_complex numbers them."""
    doc, _ = _read(structure)
    return _complex.find(doc.model())


def cdrs(structure):
    """CDRs of the TCR chains: list of dicts (chain, type, cdr, residues, length, sequence)."""
    from . import imm_cdr
    _require_anarci()
    doc, _ = _read(structure)
    out = []
    for cid, kind, k, rs in imm_cdr.cdrs(doc.model()):
        out.append(dict(chain=cid, type=kind, cdr=f"CDR{k}",
                        residues=[_io.resid(r) for r in rs], length=len(rs),
                        sequence=_chains.sequence(rs)))
    return out


def check(structure, chains=None):
    """imm_check as a list of dicts (one per check)."""
    from . import imm_check
    doc, name = _read(structure)
    return list(imm_check.rows(doc.model(), _chainset(chains), name))


def renumber(structure, out=None, tcr=True, mhc=True, peptide=True, constant=False,
             chains=None, quiet=True):
    """IMGT-number TCR, MHC and peptide chains (imm_tcrnum + imm_mhcnum + imm_pepnum).
    Returns the new file content; writes it to `out` if given."""
    doc, _ = _read(structure)
    cs = _chainset(chains)
    new = {}
    with _maybe_quiet(quiet):
        model = doc.model()
        if mhc:
            new.update(_number.mhc_ids(model, cs))
        if tcr:
            _require_anarci()
            new.update(_number.tcr_ids(model, cs, constant=constant))
        if peptide:
            new.update(_number.pep_ids(model, cs))
    return _write(doc.renumbered(new).text(), out)


def rechain(structure, out=None, layout=None, quiet=True):
    """Rename the chains of a single complex to the standard layout (imm_rechain).
    layout: {role: chain ID} replacing defaults (roles: see immtools._registry.LAYOUT)."""
    from . import imm_rechain
    option = ",".join(f"{r}:{c}" for r, c in (layout or {}).items())
    text, _ = _source(structure)
    with _maybe_quiet(quiet):
        result = imm_rechain.run(io.StringIO(text), option)
    return _write(result, out)


def select_complex(structure, n=1, out=None, quiet=True):
    """Only the chains of complex n (1-based, as imm_complex)."""
    from . import imm_complex
    text, _ = _source(structure)
    with _maybe_quiet(quiet):
        result = imm_complex.run(io.StringIO(text), str(n))
    return _write(result, out)


def structure(structure):
    """The first model as a Biopython model (chains -> residues -> atoms), from a path, an
    open file or file content (PDB or mmCIF). Handy for your own analyses after renumber()."""
    doc, _ = _read(structure)
    return doc.model()


# ------------------------------------------------------------------ helpers

def _source(structure):
    """(file content, name) from a path, an open file or a string with the content."""
    if hasattr(structure, "read"):
        return structure.read(), Path(getattr(structure, "name", "input")).stem
    s = str(structure)
    if "\n" not in s and os.path.exists(s):
        return Path(s).read_text(), Path(s).stem
    if "\n" in s:
        return s, "input"
    raise FileNotFoundError(s)


def _read(structure):
    text, name = _source(structure)
    return _doc.read(io.StringIO(text)), name


def _chainset(chains):
    if chains is None:
        return None
    return {chains} if isinstance(chains, str) else set(chains)


def _write(text, out):
    if out is not None:
        Path(out).write_text(text)
    return text


def _require_anarci():
    if not _tcr.available():
        raise RuntimeError("ANARCI is not installed: pip install anarci-pyhmmer")


class _maybe_quiet:
    """Silence the per-chain stderr summaries (and turn sys.exit into exceptions)."""

    def __init__(self, quiet):
        self.quiet = quiet
        self.buf = io.StringIO()

    def __enter__(self):
        if self.quiet:
            self.r = redirect_stderr(self.buf)
            self.r.__enter__()
        return self

    def __exit__(self, typ, val, tb):
        if self.quiet:
            self.r.__exit__(typ, val, tb)
        if typ is SystemExit:
            msg = self.buf.getvalue().strip() or "failed"
            raise ValueError(msg.replace("ERROR!! ", "")) from None
        return False
