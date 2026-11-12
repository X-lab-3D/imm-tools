# SPDX-License-Identifier: Apache-2.0
"""
Numbers the peptide 1..n (insertion codes removed). Peptides are chains of at most 30
residues without an MHC or TCR domain; with -<chain id(s)> those chains are numbered
whatever their type.

Usage:
    imm_pepnum [-<chain id(s)>] <pdb or mmcif file>

Example:
    imm_pepnum 1ao7.pdb > out.pdb
    imm_pepnum -C 1ao7.pdb
    imm_tcrnum 1ao7.pdb | imm_mhcnum | imm_pepnum > 1ao7_imgt.pdb

This program is part of imm-tools, made in the style of pdb-tools.
"""
import sys

from immtools import _io, _doc, _number


def run(fhandle, option):
    doc = _doc.read(fhandle)
    model = doc.model()
    new = _number.pep_ids(model, _io.chain_list(option))
    return doc.renumbered(new).text()


def main():
    option, fh = _io.check_input(sys.argv[1:], __doc__)
    _io.chain_list(option)
    _io.write(run(fh, option))
    fh.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
