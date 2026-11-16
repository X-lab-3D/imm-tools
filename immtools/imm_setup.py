# SPDX-License-Identifier: Apache-2.0
"""
Checks that the ANARCI backend used for TCR numbering (imm_tcrnum, imm_cdr, imm_seq,
imm_annotate, imm_check, imm_complex, imm_rechain) is installed and works. ANARCI comes
with imm-tools as the package anarci-pyhmmer (IMGT HMMs included), so nothing needs to be
set up; this tool only diagnoses problems.

Usage:
    imm_setup

This program is part of imm-tools, made in the style of pdb-tools.
"""
import sys

from immtools import _tcr


def main():
    if len(sys.argv) > 1:
        sys.stderr.write(__doc__)
        sys.exit(0 if sys.argv[1] in ("-h", "--help") else 1)
    ok, msg = _tcr.status()
    sys.stderr.write(msg + "\n")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
