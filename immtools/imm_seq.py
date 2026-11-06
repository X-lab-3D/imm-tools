# SPDX-License-Identifier: Apache-2.0
"""
Prints the sequence of every protein chain in FASTA format, with its type:
TCR_alpha / TCR_beta / TCR_gamma / TCR_delta, MHC_I, MHC_II_alpha, MHC_II_beta, B2M,
peptide, antibody_heavy / _kappa / _lambda, or other.

Usage:
    imm_seq [-<type>] <pdb or mmcif file>

    type: tcr, mhc (MHC chains and beta2m), pep, ab (antibody), or chain ids (-A,B);
    default: all chains.

Example:
    imm_seq 1ao7.pdb              # all chains
    imm_seq -mhc 1ao7.pdb         # MHC heavy chain and beta2m only
    imm_seq -tcr 1ao7.pdb > tcr.fasta
    pdb_selmodel -1 model.pdb | imm_seq -pep

Header: >name:chain|type|details. TCR identification needs the ANARCI backend
(pip install anarci-pyhmmer); without it TCR chains are reported as 'other'.

This program is part of imm-tools, made in the style of pdb-tools.
"""
import os
import sys

from immtools import _io, _doc, _chains
from immtools._registry import GROUPS as _GROUPS, family

GROUPS = {g: (lambda k, fams=fams: family(k) in fams) for g, fams in _GROUPS.items()}


def run(fhandle, option, name="input"):
    model = _doc.read(fhandle).model()
    sel = GROUPS.get(option.lower()) if option else None
    chains = None if (not option or sel) else _io.chain_list(option)
    for chain in model:
        if chains and chain.id not in chains:
            continue
        c = _chains.identify(chain)
        if c["kind"] == "none" or (sel and not sel(c["kind"])):
            continue
        desc = _chains.describe(c)
        detail = desc[len(c["kind"]):].strip(" ()")
        yield f">{name}:{chain.id}|{c['kind']}" + (f"|{detail}" if detail else "") + "\n"
        for i in range(0, len(c["seq"]), 60):
            yield c["seq"][i:i + 60] + "\n"


def main():
    option, fh = _io.check_input(sys.argv[1:], __doc__)
    if option and option.lower() not in GROUPS:
        _io.chain_list(option)
    name = os.path.splitext(os.path.basename(fh.name))[0] if fh is not sys.stdin else "stdin"
    _io.write(list(run(fh, option, name)))
    fh.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
