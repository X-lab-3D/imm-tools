# SPDX-License-Identifier: Apache-2.0
"""New residue ids for imm_tcrnum, imm_mhcnum and imm_pepnum."""
from . import _io, _chains, _tcr
from ._registry import family


def _report(chain, c, done, extra=""):
    _chains.note(f"chain {chain}: {_chains.describe(c)}{extra}" + ("" if done else " - left unchanged"))


def mhc_ids(model, chains=None, quiet=False):
    new = {}
    for ch in model:
        if chains and ch.id not in chains:
            continue
        c = _chains.identify(ch, use_tcr=False)
        if not c["mhc"]:
            if not quiet and c["kind"] != "none":
                _report(ch.id, c, False)
            continue
        if not quiet:
            _report(ch.id, c, True)
        for r, l in zip(c["residues"], c["mhc"][0]):
            new[_io.residue_key(r)] = l
    return new


def tcr_ids(model, chains=None, constant=False, quiet=False):
    new = {}
    for ch in model:
        if chains and ch.id not in chains:
            continue
        c = _chains.identify(ch)
        if not c["tcr"] or family(c["kind"]) != "tcr":     # antibodies: not (yet) here
            if not quiet and c["kind"] != "none":
                _report(ch.id, c, False)
            continue
        labels, info = _tcr.full_labels(c["seq"], constant)
        extra = ""
        if constant:
            extra = (f", constant domain {info['constant'][0]} {info['constant'][1]:.0f} bits"
                     if info["constant"] else ", no constant domain found")
        if not quiet:
            _report(ch.id, c, True, extra)
        for r, l in zip(c["residues"], labels):
            new[_io.residue_key(r)] = l
    return new


def pep_ids(model, chains=None, quiet=False):
    """Peptide chains (or the given chains) numbered 1..n."""
    new = {}
    for ch in model:
        if chains and ch.id not in chains:
            continue
        rs = _chains.protein_residues(ch)
        if not chains and not (rs and len(rs) <= _chains.PEPTIDE_MAX):
            continue
        if not chains:
            c = _chains.identify(ch)
            if c["kind"] != "peptide":
                continue
        if not quiet:
            _chains.note(f"chain {ch.id}: peptide, {len(rs)} residues numbered 1-{len(rs)}")
        for k, r in enumerate(rs, 1):
            new[_io.residue_key(r)] = (k, "")
    return new
