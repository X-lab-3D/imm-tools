<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/logo-lockup-dark.png">
    <img src="docs/images/logo-lockup.png" alt="imm-tools" width="380">
  </picture>
</p>

Consistent IMGT numbering and chain IDs for TCR–pMHC structures (TCR, MHC of any species,
β2m, peptide), so that the same residue number means the same place in every structure.

Command-line tools in the style of [pdb-tools](https://github.com/haddocking/pdb-tools)
(PDB or mmCIF in, stdout out, chainable with pipes), plus a Python interface.

**Scope:** TCR–pMHC structures (TCR α/β/γ/δ, MHC class I and II, β2m, peptide). Antibody
chains are recognised (`imm_seq`, `imm_annotate`) but not numbered or paired; the tools are
built so that antibody support (V and constant-domain numbering, H/L pairing, antigen
handling) can be added – contributions welcome; [CONTRIBUTING.md](CONTRIBUTING.md) explains
how, including a step-by-step guide for antibodies.

![The alpha1 helix of MHC class I from five species, IMGT-numbered](docs/images/alpha1_species.png)

*`imm_mhcnum` on MHC class I of five species: the same IMGT position is the same place in the
structure (chicken and carp have an insertion, 54A; position 65 is residue 63–65 of the
sequences). Bars: how many of 10 TCR–pMHC complexes touch each position
([examples/contact_map.py](examples/contact_map.py)); 65, 66, 69, 72 and 73 are the classic
TCR contacts. Made with [docs/make_figures.py](docs/make_figures.py).*

| tool | what it does |
|---|---|
| `imm_annotate` | each chain's type, complex, species, closest gene / allele and domain ranges (TSV) |
| `imm_seq` | chain sequences (FASTA) with their type |
| `imm_complex` | pick one TCR–pMHC complex out of a multi-copy file |
| `imm_rechain` | rename chains to the standard layout (A, B, M, N, 2, P), configurable |
| `imm_tcrnum` | IMGT-number the TCR (V domain; `-const` also the constant domain) |
| `imm_mhcnum` | IMGT-number the MHC chains and β2m (G and C-like domains, fixed domain offsets), any species |
| `imm_pepnum` | number the peptide 1…n |
| `imm_split` | write TCR, pMHC, MHC and peptide as separate files |
| `imm_cdr` | print or select the CDRs |
| `imm_check` | sanity-check the numbering (conserved Cys/Trp/Phe anchors, MHC numbering confidence, chain breaks) |
| `imm_setup` | check that the ANARCI backend used for TCR chains works |

Chain types: TCR α / β / γ / δ, MHC class I heavy chain, MHC class II α / β, β2m, peptide
(protein chain of ≤ 30 residues without a domain), antibody heavy / κ / λ, other.

## Install

```bash
pip install git+https://github.com/X-lab-3D/imm-tools.git
```

or from a clone:

```bash
git clone https://github.com/X-lab-3D/imm-tools.git
cd imm-tools
pip install -e .
```

This installs everything, including pdb-tools and ANARCI. If pip refuses with
"externally-managed-environment" (e.g. Homebrew Python), install into a conda or virtual
environment, or use `pipx install --include-deps git+https://github.com/X-lab-3D/imm-tools.git`.

## Typical use

```bash
imm_annotate 3qdg.pdb                         # what is in the file?

# one complex, standard chains, IMGT numbering everywhere, then check
imm_complex -1 3qdg.pdb | imm_rechain \
  | imm_tcrnum -const | imm_mhcnum | imm_pepnum > 3qdg_std_imgt.pdb
imm_check 3qdg_std_imgt.pdb

imm_split 3qdg_std_imgt.pdb                   # 3qdg_std_imgt_TCR.pdb, _pMHC, _MHC, _peptide
imm_cdr 3qdg_std_imgt.pdb                     # CDR table
imm_cdr -cdr3 3qdg_std_imgt.pdb > cdr3.pdb
imm_seq -tcr 3qdg.pdb > tcr.fasta

imm_tcrnum -const 7xyz.cif | imm_mhcnum > 7xyz_imgt.cif   # mmCIF in, mmCIF out
```

What the tables look like (`imm_annotate 3qdg.pdb`, shortened):

```
chain  type       complex  length  species       gene_allele                                  identity      domains
A      MHC_I      1        275     Homo sapiens  HLA-A*02:01 / HLA-A*02:09 / ...              100.0         G_ALPHA1:1-90;G_ALPHA2:91-182;C_ALPHA3:183-275
B      B2M        1        100     Homo sapiens  Beta-2-microglobulin [UniProt P61769] / ...  99.0          C_B2M:1-99
C      peptide    1        10
D      TCR_alpha  1        199     Homo sapiens  TRAV12-2*01 / TRAJ23*01                      97.8 / 100.0  V:1-109;C:110-199
E      TCR_beta   1        242     Homo sapiens  TRBV6-4*01 / TRBJ1-1*01                      98.9 / 100.0  V:4-115;C:116-245
```

`imm_check 3qdg_imgt.pdb` (shortened):

```
chain  type       check       status  detail
A      MHC_I      numbering   OK      275/275 residues as IMGT (100.0%)
A      MHC_I      anchors     OK      C1011 C1074 C2023 W2041 C2104
A      MHC_I      confidence  OK      lowest domain score C_ALPHA3 148 bits
D      TCR_alpha  anchors     OK      C23 W41 C104 F118
E      TCR_beta   breaks      OK      no chain breaks
```

The real output is tab-separated, starts with a `file` column, and opens in any spreadsheet.
More worked examples – a whole folder at once, a CDR3 table, a contact map by IMGT position,
PyMOL colouring – are in [examples/](examples/).

Each tool prints its help with `-h`. Options follow pdb-tools: one option, written
`-<value>`:

| tool | option |
|---|---|
| `imm_annotate`, `imm_mhcnum`, `imm_pepnum`, `imm_check` | `-A,B` chains (default: all) |
| `imm_seq` | `-tcr`, `-mhc`, `-pep` or chains |
| `imm_complex` | `-2` complex number (default 1); all complexes are listed on stderr |
| `imm_rechain` | `-tcra:D,tcrb:E` (roles tcra, tcrb, mhc, mhcb, b2m, pep) or six IDs `-A,B,M,N,2,P` |
| `imm_tcrnum` | `-const`, chains, or `-const:D,E` |
| `imm_split` | parts `-tcr,pmhc,mhc,pep` (default: all) |
| `imm_cdr` | none (table) or `-cdr1`, `-cdr3`, `-a3,b3`, … (PDB selection) |

What the tools change: the numbering tools change only the residue number and insertion code
columns of ATOM / HETATM / ANISOU / TER records, `imm_rechain` only the chain column; the
residue references in HELIX, SHEET, SSBOND, LINK and CISPEP records follow, so viewers such as
PyMOL that read them still show the right secondary structure. All other lines and columns
are kept (SEQRES and other header records are not updated). Waters and ligands keep their
numbers unless they collide with a new residue number (then they move after the chain's last
residue). With several models, the first model is
analysed and the result applied to all. Files with bare `TER` lines break some pdb-tools;
`pdb_tidy` fixes them.

**PDB and mmCIF.** Every tool reads both (the format is recognised from the content) and
writes the same format as its input. mmCIF is needed for large entries (more than 99,999
atoms, or chain IDs longer than one character), which have no PDB-format file. In mmCIF
files the author fields change (`auth_asym_id`, `auth_seq_id`, `pdbx_PDB_ins_code`, and the
same fields in `struct_conf`, `struct_sheet_range`, `struct_conn`, `struct_mon_prot_cis`,
`pdbx_poly_seq_scheme`); `label_*` fields and all other categories are kept. The pdb-tools
commands (`pdb_*`) are PDB-only.

## Python

```python
import immtools

n = immtools.number_mhc("GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEY...")
n.chain_type                  # 'MHC_I'
n.labels[:3]                  # ['1', '2', '3']
n.domains                     # [('G_ALPHA1', 0, 88, 150.7), ('G_ALPHA2', 90, 175, 142.8), ...]
immtools.number_tcr(seq, constant=True)   # .genes: species, V / J genes

immtools.annotate("3qdg.pdb")              # list of dicts, as imm_annotate
immtools.complexes("1d9k.pdb")             # [{'mhc': 'C', 'mhcb': 'D', 'pep': 'P', ...}, ...]
immtools.cdrs("3qdg.pdb")                  # chain, cdr, residues, sequence
immtools.check("3qdg_imgt.pdb")            # list of dicts, as imm_check

text = immtools.renumber("3qdg.cif", out="3qdg_imgt.cif", constant=True)
text = immtools.rechain(immtools.select_complex("1d9k.pdb", 2), out="1d9k_2.pdb")

model = immtools.structure(text)           # Biopython model, for your own analyses
model["D"][104].get_resname()              # 'CYS': the conserved Cys at IMGT 104
```

Structures can be a file path (PDB or mmCIF), an open file or the file content as a
string; functions that change a structure return the new content and write it to `out` if
given. These functions (`help(immtools)`) are the stable interface; modules whose names start
with an underscore are internal.

## Complexes and chain layout

`imm_complex`, `imm_rechain`, `imm_split` and `imm_annotate` group chains into complexes by
type and contacts (CA–CA ≤ 8 Å): MHC heavy chain (or class II α) → β2m (or class II β) → the
peptide touching the groove → the TCR α/β pair touching the pMHC. TCRs without pMHC and lone
MHC chains form complexes of their own.

Standard layout (`imm_rechain`): A TCR α (or γ), B TCR β (or δ), M MHC class I heavy chain /
class II α, N MHC class II β, 2 β2m, P peptide. β2m gets a digit so that H and L stay free for
antibody heavy and light chains, and because lower-case IDs are matched case-insensitively by
some viewers (PyMOL selections by default).

## Numbering

**TCR** (`imm_tcrnum`): IMGT numbering of the V domain by ANARCI. Without `-const`, residues after
the V domain continue after the highest IMGT number. With `-const`, the constant domain gets
IMGT C-LIKE numbering + 1000 (1001–1128; TRAC, TRBC, TRGC, TRDC; profile HMMs built from
IMGT/DomainDisplay, all species), residues after it 2001, 2002, …. Residues before the V
domain: 0, −1, … . A second V domain in the same chain (single-chain constructs) gets +1000.
Note: IMGT/3Dstructure-DB files number the TRAC constant domain 2001–2128; here every constant
domain is 1001–1128.

**MHC** (`imm_mhcnum`): IMGT unique numbering (G-DOMAIN 1–92, C-LIKE 1–128, IMGT insertion
codes such as 49A, 61A, 72A, 84A/85A) plus 1000 × the domain's place in the chain:

| chain | domain 1 | domain 2 | domain 3 |
|---|---|---|---|
| class I heavy chain | α1 1–92 | α2 1001–1092 | α3 2001–2128 |
| β2-microglobulin | 1–128 | | |
| class II α | α1 1–92 | α2 1001–1128 | |
| class II β | β1 1–92 | β2 1001–1128 | |

Residues before the first domain: 0, −1, …; after the last one: the next thousand. Chain
breaks are found from the coordinates, so missing loops do not shift the numbering.

MHC numbering uses two sets of profile HMMs with IMGT labels: a **broad** set covering all
jawed-vertebrate lineages (fish, amphibians, reptiles, birds, mammals) and a **mammal-focused**
set. Every domain is searched with both and numbered by the set that matches it best (higher
HMMER bit score); `imm_check` lists the residues the two sets number differently, which are
the least certain ones. The models (`immtools/models/mhc_broad`, `mhc_mammal`) were trained on
MHC sequences from IPD-IMGT/HLA, IPD-MHC and UniProt, with IMGT structures and AlphaFold DB
models as references; a description of the method will follow in a publication. The model version is in
`immtools/models/MODELS_VERSION.json`.

**Peptide** (`imm_pepnum`): 1…n.

## Annotation

`imm_annotate` columns: file, chain, type, complex, length, species, gene_allele,
identity, domains.

* TCR: species and V / J germline genes from ANARCI, with the V and J identity.
* MHC: the closest allele / gene by the peptide-binding region (class I α1+α2, class II α1
  or β1; β2m: whole domain) among ~57,000 distinct regions from IPD-IMGT/HLA, IPD-MHC and
  UniProt, with % identity over that region. Alleles that share the region are listed
  together ("+N" more), so the 2-field allele is only resolved where the peptide-binding
  region differs.
* domains: `name:first-last` residue, in the numbering of the input file.

## How accurate is it? Please help test

Share of residues numbered as in the reference (IMGT/3Dstructure-DB unless noted):

| chains | test set | agreement |
|---|---|---|
| MHC class I | human and mouse structures¹ | 100 % |
| β2m | human and mouse structures¹ | 100 % |
| MHC class II α | human and mouse structures¹ | 99.96 % |
| MHC class II β | human and mouse structures¹ | 99.74 % |
| MHC, other species (chicken, macaque, pig, …) | structures, species held out | 98.5–100 % |
| MHC class I, grass carp | structures, species held out | 93.6 % |
| MHC, birds, reptiles, amphibians | AlphaFold models², lineage held out | 94–99 % |
| MHC class I / class II β, fish | AlphaFold models², lineage held out | 92.4 / 93.6 % |
| TCR constant domain, TRAC / TRBC | 143 distinct TCR chains, mostly human and mouse | 98.8 / 98.6 % |

¹ These structures are in the training data. ² No IMGT numbering exists for these; the
reference is the structurally equivalent position.

* Fish chains are the weakest; errors concentrate in loops and at domain ends. γδ constant
  domains are untested.
* Where IMGT numbers the same loop length differently in different species, one convention
  is used; where the structure clearly disagrees with an IMGT label, the structure wins. All
  such residues are listed in `immtools/models/mhc_broad/differences_from_imgt.tsv`.
* Not specifically tested: non-classical MHC (MR1, CD1, HFE, MIC, …), single-chain TCRs,
  peptides tethered to the MHC chain.

Please report structures where the output looks wrong (PDB ID, chain, residue), with the
stderr summary and the `imm_check` table.

## Tests

```bash
pip install pytest && pytest tests/        # unit tests + end-to-end tests of every tool
```

## Acknowledgements and third-party components

imm-tools builds on the following software and data. Please cite them when you use the
corresponding tools.

Software
* **ANARCI** – TCR (and antibody) V-domain numbering and chain-type / germline assignment
  (`imm_tcrnum`, `imm_cdr`, `imm_seq`, `imm_annotate`, `imm_check`). J. Dunbar and C. M. Deane,
  *Bioinformatics* 32, 298–300 (2016). © Oxford Protein Informatics Group, BSD-3-Clause.
  imm-tools uses the **pyhmmer port of ANARCI by David Prihoda**
  (https://github.com/prihoda/ANARCI), packaged with the IMGT HMMs as `anarci-pyhmmer` from
  the fork https://github.com/X-lab-3D/ANARCI-pyhmmer-fork. The fork includes MUSCLE (public
  domain), used only to rebuild the HMMs.
* **pyhmmer** – M. Larralde and G. Zeller, *Bioinformatics* 39, btad214 (2023); MIT. Python
  bindings to **HMMER** (S. R. Eddy, *PLoS Comput. Biol.* 7, e1002195 (2011)).
* **Biopython** – P. J. A. Cock et al., *Bioinformatics* 25, 1422–1423 (2009).
* **gemmi** – M. Wojdyr, *J. Open Source Softw.* 7, 4200 (2022); MPL-2.0. mmCIF reading and
  writing.
* **pdb-tools** – J. P. G. L. M. Rodrigues et al., *F1000Research* 7, 1961 (2018); the design
  of these tools follows pdb-tools.
* **TM-align** via **tmtools** – Y. Zhang and J. Skolnick, *Nucleic Acids Res.* 33, 2302–2309
  (2005); used only to build the MHC models, not at run time.

Data used to build the models and reference sequences
* **IMGT®**, the international ImMunoGeneTics information system® (https://www.imgt.org):
  IMGT/3Dstructure-DB and IMGT/DomainDisplay – IMGT unique numbering, MHC and TCR domain
  alignments. M.-P. Lefranc et al.; F. Ehrenmann, Q. Kaas and M.-P. Lefranc, *Nucleic Acids
  Res.* 38, D301–D307 (2010).
* **IPD-IMGT/HLA** (https://www.ebi.ac.uk/ipd/imgt/hla/) – J. Robinson et al.
* **IPD-MHC** (https://www.ebi.ac.uk/ipd/mhc/) – G. Maccari et al.
* **UniProt** (https://www.uniprot.org) – The UniProt Consortium; InterPro domain annotations.
* **AlphaFold Protein Structure Database** (https://alphafold.ebi.ac.uk) – M. Varadi et al.,
  *Nucleic Acids Res.* (2022, 2024); models under CC-BY-4.0.

IMGT® is a registered trademark of CNRS. imm-tools is not affiliated with or endorsed by
IMGT, IPD, UniProt or EMBL-EBI.

## Licence

Apache-2.0 (see LICENSE and NOTICE). TCR V-domain numbering uses ANARCI (BSD-3-Clause). The
MHC and TCR-constant models and the MHC reference sequences are derived from
IMGT/3Dstructure-DB, IMGT/DomainDisplay, IPD-IMGT/HLA, IPD-MHC, UniProt and AlphaFold DB data;
check their terms before redistributing them. The scripts that built the TCR constant-domain
models and the MHC reference sequences are in `dev/`.
