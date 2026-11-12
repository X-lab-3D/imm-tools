# SPDX-License-Identifier: Apache-2.0
"""imm-tools: consistent IMGT numbering and chain IDs for TCR-pMHC structures (pdb-tools style).

    imm_annotate  chain type, complex, species, closest gene / allele, domain ranges (TSV)
    imm_seq       chain sequences (FASTA) with their type
    imm_complex   pick one TCR-pMHC complex out of a multi-copy file
    imm_rechain   rename chains to the standard layout (A, B, M, N, 2, P), configurable
    imm_tcrnum    IMGT-number TCR chains (V domain; -const also the constant domain)
    imm_mhcnum    IMGT-number MHC chains and beta2m (any species)
    imm_pepnum    number the peptide 1..n
    imm_split     write TCR, pMHC, MHC and peptide as separate files
    imm_cdr       print or select CDRs
    imm_check     sanity-check the numbering
    imm_setup     check the ANARCI backend (anarci-pyhmmer)
"""
__version__ = "0.4.0"

from .api import (Numbering, number_mhc, number_tcr, annotate, complexes, cdrs, check,  # noqa
                  renumber, rechain, select_complex, structure)
