# SPDX-License-Identifier: Apache-2.0
"""
IMGT-numbers the TCR chains (alpha, beta, gamma, delta) of a PDB or mmCIF file: the variable domain
with ANARCI and, with -const, the constant domain too. Other chains and all other records
are left unchanged.

Usage:
    imm_tcrnum [-const | -<chain id(s)> | -const:<chain id(s)>] <pdb or mmcif file>

Example:
    imm_tcrnum 1ao7.pdb > 1ao7_tcr.pdb
    imm_tcrnum -const 1ao7.pdb                   # V and C domain
    imm_tcrnum -const:D,E 1ao7.pdb               # V and C domain of chains D and E

Numbering:
    V domain   IMGT 1-128 (CDR3 insertions 111A.., 112A..); a second V domain in the same
               chain (single-chain constructs) +1000
    default    residues after the V domain continue after the highest IMGT number;
    -const     constant domain IMGT C-LIKE 1001-1128 (TRAC, TRBC, TRGC, TRDC; profile HMMs
               from IMGT/DomainDisplay, all species); residues between V and C continue
               the V domain with letters; residues after C: 2001, 2002, ...
    residues before the V domain: 0, -1, ...
Note: IMGT/3Dstructure-DB files number the TRAC constant domain 2001-2128; here every
constant domain is 1001-1128.

Needs ANARCI (anarci-pyhmmer, installed with imm-tools; imm_setup checks it).

This program is part of imm-tools, made in the style of pdb-tools.
"""
import sys

from immtools import _io, _doc, _number, _tcr


def parse(option):
    const = option == "const" or option.startswith("const:")
    chains = option.split(":", 1)[1] if option.startswith("const:") else ("" if const else option)
    return const, _io.chain_list(chains)


def run(fhandle, option):
    const, chains = parse(option)
    doc = _doc.read(fhandle)
    new = _number.tcr_ids(doc.model(), chains, constant=const)
    return doc.renumbered(new).text()


def main():
    option, fh = _io.check_input(sys.argv[1:], __doc__)
    parse(option)
    _tcr.require()
    _io.write(run(fh, option))
    fh.close()
    sys.exit(0)


if __name__ == "__main__":
    main()
