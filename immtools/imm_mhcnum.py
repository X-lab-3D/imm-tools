# SPDX-License-Identifier: Apache-2.0
"""
IMGT-numbers the MHC chains (class I heavy chain, class II alpha and beta) and beta2m of a
PDB or mmCIF file. Other chains and all other records are left unchanged.

Usage:
    imm_mhcnum [-<chain id(s)>] <pdb or mmcif file>

Example:
    imm_mhcnum 1ao7.pdb > 1ao7_mhc.pdb
    imm_mhcnum -A,B 1ao7.pdb                     # only chains A and B
    imm_tcrnum 1ao7.pdb | imm_mhcnum | imm_pepnum > 1ao7_imgt.pdb

Numbering: IMGT position per domain (G-DOMAIN 1-92, C-LIKE 1-128, IMGT insertion codes)
plus 1000 x the domain's place in the chain:
    class I heavy chain   alpha1 1-92, alpha2 1001-1092, alpha3 2001-2128
    class II alpha        alpha1 1-92, alpha2 1001-1128
    class II beta         beta1 1-92, beta2 1001-1128
    beta2m                1-128
Residues before the first domain: 0, -1, ...; after the last domain: the next thousand.
Any species: two sets of profile HMMs (all jawed-vertebrate lineages, and mammal-focused);
each domain is numbered by the set that matches it best. Chain breaks are detected from the coordinates.

Waters and ligands keep their numbers unless they collide with a new residue number; then
they move after the chain's last residue. Residue references in HELIX, SHEET, SSBOND, LINK
and CISPEP records are renumbered too, so viewers show the right secondary structure. A per-chain summary is printed to stderr.

This program is part of imm-tools, made in the style of pdb-tools.
"""
import sys

from immtools import _io, _doc, _number


def run(fhandle, option):
    doc = _doc.read(fhandle)
    model = doc.model()
    new = _number.mhc_ids(model, _io.chain_list(option))
    return doc.renumbered(new).text()


def main():
    option, fh = _io.check_input(sys.argv[1:], __doc__)
    _io.chain_list(option)
    _io.write(run(fh, option))
    fh.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
