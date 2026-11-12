#!/usr/bin/env python3
"""
Score docked TCR-pMHC models against the bound (reference) structure.

Docking starts from unbound structures, so models and reference usually differ in chain IDs,
residue numbering, missing residues and even constructs. This script standardises both with
imm-tools (first complex, IMGT numbering, standard chain IDs), so residues are matched simply
by (chain, number, insertion code), and then reports

  L-RMSD   TCR V domains (backbone), after superposing on the pMHC (MHC + peptide)
  I-RMSD   interface residues (backbone, any heavy atom within 10 A of the partner in the
           reference), after superposing on them
  fnat     fraction of the reference TCR-pMHC residue contacts (heavy atoms within 5 A)
           present in the model
  DockQ    (fnat + 1/(1+(I_RMSD/1.5)^2) + 1/(1+(L_RMSD/8.5)^2)) / 3  (Basu & Wallner 2016)

Only residues present in both structures are used. The numbers agree with the DockQ program
(v2) run on the same standardised files with the TCR and the pMHC each merged into one chain.

Usage:
    python docking_rmsd.py reference.pdb model1.pdb [model2.pdb ...]
"""
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import NeighborSearch, Superimposer
from Bio.PDB.Polypeptide import is_aa

import immtools

TCR, PMHC = ("A", "B"), ("M", "P")         # standard chain IDs (class I; class II adds N)
BACKBONE = ("N", "CA", "C", "O")


def standardise(path):
    """First complex, IMGT-numbered (TCR V + C), standard chain IDs -> Biopython model."""
    text = immtools.select_complex(path, 1)
    text = immtools.renumber(text, constant=True)
    return immtools.structure(immtools.rechain(text))


def residues(model, chains):
    out = {}
    for c in chains:
        if c not in model:
            continue
        for r in model[c]:
            if not is_aa(r):
                continue
            if c in TCR and r.id[1] >= 1000:          # TCR: variable domains only
                continue
            out[(c, r.id[1], r.id[2])] = r
    return out


def atoms(res, keys):
    """Backbone atoms present in both, as paired lists."""
    a, b = [], []
    for k in keys:
        for name in BACKBONE:
            if name in res[0][k] and name in res[1][k]:
                a.append(res[0][k][name])
                b.append(res[1][k][name])
    return a, b


def rmsd_after_fit(fit, measure):
    """Superpose model onto reference by `fit` pairs, RMSD over `measure` pairs."""
    sup = Superimposer()
    sup.set_atoms(fit[0], fit[1])                     # fixed = reference, moving = model
    rot, tran = sup.rotran
    ref = np.array([a.coord for a in measure[0]])
    mod = np.array([a.coord for a in measure[1]]) @ rot + tran
    return float(np.sqrt(((ref - mod) ** 2).sum(1).mean()))


def contacts(model, cutoff):
    heavy = [a for c in TCR + PMHC if c in model for r in model[c] if is_aa(r)
             for a in r if a.element != "H"]
    pairs = set()
    for a, b in NeighborSearch(heavy).search_all(cutoff):
        ra, rb = a.get_parent(), b.get_parent()
        ca, cb = ra.get_parent().id, rb.get_parent().id
        if (ca in TCR) == (cb in TCR):
            continue
        key = lambda r: (r.get_parent().id, r.id[1], r.id[2])
        pairs.add(tuple(sorted((key(ra), key(rb)))))
    return pairs


def score(ref, mod):
    rt, mt = residues(ref, TCR), residues(mod, TCR)
    rp, mp = residues(ref, PMHC), residues(mod, PMHC)
    both = lambda x, y: [k for k in x if k in y]
    tcr_keys, pmhc_keys = both(rt, mt), both(rp, mp)
    allr, allm = {**rt, **rp}, {**mt, **mp}

    lrmsd = rmsd_after_fit(atoms((allr, allm), pmhc_keys), atoms((allr, allm), tcr_keys))

    ref_contacts10 = contacts(ref, 10.0)
    iface = sorted({k for pair in ref_contacts10 for k in pair} & set(allm))
    irmsd = rmsd_after_fit(atoms((allr, allm), iface), atoms((allr, allm), iface))

    native = contacts(ref, 5.0)
    native_shared = {p for p in native if p[0] in allm and p[1] in allm}
    fnat = len(native_shared & contacts(mod, 5.0)) / max(1, len(native_shared))
    dockq = (fnat + 1 / (1 + (irmsd / 1.5) ** 2) + 1 / (1 + (lrmsd / 8.5) ** 2)) / 3
    return dict(L_RMSD=lrmsd, I_RMSD=irmsd, fnat=fnat, DockQ=dockq, tcr_residues=len(tcr_keys),
                pmhc_residues=len(pmhc_keys), interface_residues=len(iface))


def main(reference, models):
    ref = standardise(reference)
    print("model\tL_RMSD\tI_RMSD\tfnat\tDockQ\ttcr_residues\tpmhc_residues\tinterface_residues")
    for m in models:
        s = score(ref, standardise(m))
        print(f"{Path(m).name}\t{s['L_RMSD']:.2f}\t{s['I_RMSD']:.2f}\t{s['fnat']:.2f}\t{s['DockQ']:.3f}\t"
              f"{s['tcr_residues']}\t{s['pmhc_residues']}\t{s['interface_residues']}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2:])
