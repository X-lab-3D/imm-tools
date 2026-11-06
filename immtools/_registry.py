# SPDX-License-Identifier: Apache-2.0
"""All chain types, roles, layouts and per-type rules in one place.

To support a new kind of molecule (e.g. antibodies), extend the tables here first; the
tools read them. See CONTRIBUTING.md, "Adding antibodies".
"""

# --------------------------------------------------------------------------- chain types
# kind -> family. Families: "tcr", "antibody", "mhc", "peptide", "other".
KINDS = {
    "TCR_alpha": "tcr", "TCR_beta": "tcr", "TCR_gamma": "tcr", "TCR_delta": "tcr",
    "antibody_heavy": "antibody", "antibody_kappa": "antibody", "antibody_lambda": "antibody",
    "MHC_I": "mhc", "MHC_II_alpha": "mhc", "MHC_II_beta": "mhc", "B2M": "mhc",
    "peptide": "peptide", "other": "other",
}

# ANARCI chain type -> kind
ANARCI_KIND = {"A": "TCR_alpha", "B": "TCR_beta", "G": "TCR_gamma", "D": "TCR_delta",
               "H": "antibody_heavy", "K": "antibody_kappa", "L": "antibody_lambda"}

# MHC HMM chain type -> kind
MHC_KIND = {"MHC_I_ALPHA": "MHC_I", "MHC_II_ALPHA": "MHC_II_alpha",
            "MHC_II_BETA": "MHC_II_beta", "B2M": "B2M"}

PEPTIDE_MAX = 30   # residues; shorter protein chains without a domain hit are peptides

# imm_seq -<group>
GROUPS = {"tcr": ("tcr",), "mhc": ("mhc",), "pep": ("peptide",), "ab": ("antibody",)}

# --------------------------------------------------------------------------- roles
# Role of a chain inside a complex. Reserved for antibodies (not yet assigned by any
# assembler): abh (heavy), abl (light), ag (antigen).
ROLES = ("tcra", "tcrb", "mhc", "mhcb", "b2m", "pep", "abh", "abl", "ag")
ROLE_OF_KIND = {"MHC_I": "mhc", "MHC_II_alpha": "mhc", "MHC_II_beta": "mhcb", "B2M": "b2m",
                "peptide": "pep", "TCR_alpha": "tcra", "TCR_gamma": "tcra",
                "TCR_beta": "tcrb", "TCR_delta": "tcrb",
                "antibody_heavy": "abh", "antibody_kappa": "abl", "antibody_lambda": "abl"}

# imm_rechain default chain IDs (one character, unique). H / L kept free for antibodies.
LAYOUT = {"tcra": "A", "tcrb": "B", "mhc": "M", "mhcb": "N", "b2m": "2", "pep": "P",
          "abh": "H", "abl": "L", "ag": "G"}

# imm_split parts: option name -> (file suffix, roles)
PARTS = {"tcr": ("TCR", ("tcra", "tcrb")),
         "pmhc": ("pMHC", ("mhc", "mhcb", "b2m", "pep")),
         "mhc": ("MHC", ("mhc", "mhcb", "b2m")),
         "pep": ("peptide", ("pep",))}

# --------------------------------------------------------------------------- CDRs
# IMGT CDR ranges of V domains (IMGT positions, any receptor)
CDR = {1: (27, 38), 2: (56, 65), 3: (105, 117)}
# kind -> prefix used in imm_cdr options (a1, b3, ...)
CDR_SIDE = {"TCR_alpha": "a", "TCR_gamma": "a", "TCR_beta": "b", "TCR_delta": "b"}

# --------------------------------------------------------------------------- imm_check
# conserved residues at IMGT positions (file numbering), per kind or "V" for any V domain
ANCHORS = {
    "V": [(23, "C"), (41, "W"), (104, "C"), (118, "FW")],
    "MHC_I": [(1011, "C"), (1074, "C"), (2023, "C"), (2041, "W"), (2104, "C")],
    "MHC_II_alpha": [(1023, "C"), (1041, "W"), (1104, "C")],
    "MHC_II_beta": [(11, "C"), (74, "C"), (1023, "C"), (1041, "W"), (1104, "C")],
    "B2M": [(23, "C"), (104, "C")],
}


def family(kind):
    return KINDS.get(kind, "other")
