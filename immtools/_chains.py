# SPDX-License-Identifier: Apache-2.0
"""Chain identification: TCR / antibody (ANARCI), MHC and beta2m (MHC HMMs), peptide."""
import sys

from ._registry import ANARCI_KIND as TCR_KIND, MHC_KIND, PEPTIDE_MAX  # noqa: F401


def protein_residues(chain):
    from Bio.PDB.Polypeptide import is_aa
    return [r for r in chain if is_aa(r, standard=False)]


def sequence(residues):
    from Bio.Data.PDBData import protein_letters_3to1
    return "".join(protein_letters_3to1.get(r.get_resname(), "X") for r in residues)


def identify(chain, use_tcr=True):
    """dict(chain, kind, residues, seq, tcr=(labels, domains) | None, mhc=(labels, info) | None)

    kind: TCR_alpha/beta/gamma/delta, antibody_heavy/kappa/lambda, MHC_I, MHC_II_alpha,
    MHC_II_beta, B2M, peptide, other."""
    from . import _mhc, _tcr
    rs = protein_residues(chain)
    seq = sequence(rs)
    out = dict(chain=chain.id, kind="other", residues=rs, seq=seq, tcr=None, mhc=None)
    if not rs:
        out["kind"] = "none"
        return out
    if len(rs) >= 60:
        ca = [r for r in rs if "CA" in r]
        if len(ca) == len(rs):
            try:
                out["mhc"] = _mhc.number_residues(rs)
            except ValueError:
                pass
        else:                                  # residues without CA: sequence only
            try:
                out["mhc"] = _mhc.number_sequence(seq)
            except ValueError:
                pass
        if out["mhc"]:
            out["kind"] = MHC_KIND[out["mhc"][1]["chain_type"]]
            return out
        if use_tcr and _tcr.available():
            out["tcr"] = _tcr.number_chain(seq)
            if out["tcr"]:
                out["kind"] = TCR_KIND.get(out["tcr"][1][0]["chain_type"], "other")
    elif len(rs) <= PEPTIDE_MAX:
        out["kind"] = "peptide"
    return out


def identify_all(model, use_tcr=True):
    """{chain id: identify(chain)} in file order."""
    return {ch.id: identify(ch, use_tcr) for ch in model}


def describe(c):
    """One-line description of an identified chain."""
    if c["mhc"]:
        info = c["mhc"][1]
        doms = ", ".join(f"{m} {sc:.0f} bits" for m, s, e, sc in info["domains"])
        return f"{c['kind']} ({doms})"
    if c["tcr"]:
        d = c["tcr"][1][0]
        sc = f"{d['bitscore']:.0f} bits" if d.get("bitscore") is not None else ""
        return f"{c['kind']} ({d['species']}, {sc})"
    return c["kind"]


def note(msg):
    sys.stderr.write(msg + "\n")
