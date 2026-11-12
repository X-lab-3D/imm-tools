# Contributing to imm-tools

Contributions are welcome: bug reports with a PDB ID, new tools, and support for more
molecules (antibodies and nanobodies are the obvious next step – see below).

## Ground rules (pdb-tools style)

Every tool:

* does one job and is a single module `immtools/imm_<name>.py` with a `main()`;
* reads a PDB or mmCIF file given as argument **or stdin** (`_doc.read`), writes the same
  format to **stdout** (except `imm_split`,
  which writes files), and prints summaries / warnings to **stderr**;
* takes **at most one option**, written `-<value>` (parse it with `_io.check_input`);
  combine values inside the option if needed (`-const:D,E`);
* has its help text as the module docstring (`-h` prints it): what it does, usage, examples;
* changes only what it is about: numbering tools only the residue number and insertion code
  columns, `imm_rechain` only chain IDs. Edit through the document methods
  (`doc.renumbered`, `doc.renamed`, `doc.selected`, `doc.residues`), which work for PDB and
  mmCIF and keep secondary-structure and link records consistent (IMGT's descending
  insertion codes need care, see `_io.ss_runs`);
* exits with 0 on success, 1 with an `ERROR!!` message on bad input.

## Code map

| module | what is in it |
|---|---|
| `api.py` | the public Python interface (`import immtools`); keep it stable |
| `_registry.py` | **all chain types, roles, chain-ID layout, split parts, CDR ranges, check anchors** |
| `_chains.py` | `identify(chain)`: chain type from the MHC HMMs, ANARCI, or length (peptide) |
| `_complex.py` | grouping chains into complexes; one **assembler** function per kind of complex |
| `_number.py` | new residue ids for the numbering tools |
| `_doc.py` | PDB / mmCIF documents: read, renumber, rename, select, write |
| `_io.py` | option handling, PDB line-level editing |
| `_mhc.py` | MHC numbering engine (profile HMMs in `models/mhc_broad`, `models/mhc_mammal`) |
| `_tcr.py` | ANARCI V-domain numbering; TCR constant-domain HMMs (`models/tcr_constant`) |
| `_annot.py` | closest MHC allele (`models/mhc_refs.fasta.gz`) |
| `dev/` | scripts that built the models |
| `tests/` | pytest; PDB files in `tests/data/` |
| `examples/` | worked examples (scripts using the tools and the Python API); new ones welcome |
| `docs/` | README figure (`make_figures.py`) and images |

## Adding a tool

1. Copy a small tool (e.g. `imm_pepnum.py`) to `immtools/imm_<name>.py`; write the docstring
   first.
2. Register it in `pyproject.toml` under `[project.scripts]`, and list it in
   `immtools/__init__.py` and the README table.
3. If it is useful from Python, add a function to `api.py` (and `__all__`) that returns
   plain data (lists of dicts, strings) and a test in `tests/test_api.py`.
4. Add tests: end-to-end tests of the command in `tests/test_tools.py` (run it with
   `subprocess`, check the output), unit tests of new functions in `tests/test_units.py`
   (small synthetic inputs). GitHub Actions runs both on Linux and macOS for every push and
   pull request.
5. `pip install -e . && pytest tests/`

## Adding antibodies (and nanobodies, scFv, …)

Antibody chains are already **recognised**: ANARCI classifies V domains as TCR α/β/γ/δ or
antibody heavy / κ / λ, and `imm_seq` / `imm_annotate` report them with species and V / J
genes. They are not yet grouped into complexes, renumbered, renamed or split. The steps
below add that; each step is independent and can be a separate pull request.

**Conventions already decided** (please keep them, so the TCR and antibody sides fit):

* roles `abh` (heavy), `abl` (light), `ag` (antigen) – in `_registry.ROLES`;
* default chain IDs H (heavy), L (light), G (antigen) – in `_registry.LAYOUT`; β2m is `2`, so
  H / L do not clash with the pMHC layout (A, B, M, N, 2, P);
* the same tools work for every receptor type; only numbering has a tool per receptor
  (`imm_tcrnum`, `imm_mhcnum`, and for antibodies `imm_abnum`);
* chain types are detected automatically from sequence (no `-type` option);
* domain offsets as for TCR / MHC: V domain IMGT 1–128, constant domains + 1000 per domain.

### 1. Complexes – an assembler

Write a function that takes the `Context` and returns complexes (dicts role → chain ID), and
add it to `_complex.ASSEMBLERS`. Assemblers run in order and only see chains not used by
earlier ones. A working starting point (Fab + antigen, nanobody + antigen):

```python
def antibody_antigen(ctx):
    """Antibody (H + L, or a nanobody: H alone) + the antigen chain it touches most."""
    out = []
    for h in ctx.chains("abh"):
        ctx.use(h)
        cx = {"abh": h}
        ctx.take("abl", [h], cx)                       # light chain with most contacts
        paratope = [cx[r] for r in ("abh", "abl") if r in cx]
        cands = ctx.free_chains(("other", "peptide"))  # antigen: any other protein chain
        best = max(cands, key=lambda c: sum(ctx.n(c, p) for p in paratope), default=None)
        if best and sum(ctx.n(best, p) for p in paratope) >= MIN_CONTACTS:
            cx["ag"] = best
            ctx.use(best)
        out.append(cx)
    for light in ctx.chains("abl"):                    # light chains without a heavy chain
        ctx.use(light)
        out.append({"abl": light})
    return out

ASSEMBLERS = [tcr_pmhc, antibody_antigen, lone_b2m, tcr_only]
```

On 1VFB this gives `abh=B abl=A ag=C`; on 3OGO (four nanobody–antigen copies) one complex
per nanobody. Things to decide while doing it:

* antigens made of several chains (assign all chains above a contact threshold?);
* antibodies binding a pMHC (TCR-mimic antibodies): the pMHC is used by `tcr_pmhc` first, so
  either let `tcr_pmhc` accept an antibody instead of a TCR, or add the antibody to that
  complex afterwards;
* single-chain constructs (scFv: H and L V domains in one chain; ANARCI returns two domains
  – `_tcr.number_chain` already offsets the second one by 1000).

Once complexes exist, `imm_complex`, `imm_rechain` (H, L, G) and `imm_annotate` (complex
column) work without further changes. For `imm_split`, add parts to `_registry.PARTS`, e.g.
`"ab": ("antibody", ("abh", "abl"))` and `"ag": ("antigen", ("ag",))`.

### 2. Numbering – `imm_abnum`

Copy `imm_tcrnum.py` to `imm_abnum.py` and select antibody chains instead of TCR chains
(`_number.tcr_ids` skips chains whose `_registry.family(kind)` is not `"tcr"`; write
`ab_ids` the same way for `"antibody"`). ANARCI supports the IMGT, Kabat, Chothia, Martin
and AHo schemes (`_tcr.number_chain(seq, scheme=...)`); antibody users will want an option
such as `-kabat`. Note that Kabat / Chothia insertion codes ascend, IMGT's CDR3 ones descend –
`_io.apply_numbering` handles both.

Constant domains (CH1, CL κ, CL λ, hinge, CH2, CH3) with IMGT C-LIKE numbering: build HMMs
like the TCR ones with `dev/build_tcr_constant.py` from IMGT/DomainDisplay with receptor
type IG (add the IG domain labels to its `KIND` map), put them in `models/` and label them
like `_tcr.number_constant`. Decide and document the offsets (e.g. CH1 1001–, CH2 2001–,
CH3 3001–).

### 3. CDRs, checks, annotation

* `imm_cdr`: add `"antibody_heavy": "h", "antibody_kappa": "l", "antibody_lambda": "l"` to
  `_registry.CDR_SIDE` (options `-h3`, `-l1`, …) and filter by family in `imm_cdr.cdrs`.
  For other schemes, CDR definitions differ (Kabat, Chothia, North): add them per scheme.
* `imm_check`: the V-domain anchors (Cys23, Trp41, Cys104, Phe/Trp118) apply to antibodies in
  IMGT numbering already (`_registry.ANCHORS["V"]`); add constant-domain anchors once they
  are numbered.
* `imm_annotate`: V / J genes and species come from ANARCI already. Its germline HMMs cover
  antibodies of human, mouse, rat, rabbit, rhesus and pig (TCRs: human and mouse only), so
  e.g. camelid nanobodies are reported as the closest of those (3OGO: human IGHV3-11, 71 %).
  Camelid germlines, and the isotype (IgG1, IgA, κ / λ) from the constant domain, would be
  useful additions.

### 4. Tests

Add a few public PDB entries to `tests/data/` (a Fab–antigen, a nanobody–antigen, an scFv, a
TCR-mimic antibody–pMHC) and test: chain types, complexes, `imm_rechain` layout, numbering of
conserved positions, CDR sequences, and that `imm_tcrnum` still leaves antibody chains alone.
Keep every existing test passing.

## Pull requests

* one topic per pull request; describe what you tested (PDB IDs);
* `pytest tests/` passes;
* README updated (tool table, options, numbering scheme) and, for new data or software,
  the acknowledgements and NOTICE.
