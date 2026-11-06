"""Unit tests of the internal modules (fast, small synthetic inputs).

The end-to-end tests of the command-line tools are in test_tools.py.
Run with:  pytest tests/
"""
import io
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from immtools import _io, _registry, _complex, _mhc, _tcr, _annot          # noqa: E402
from immtools import imm_rechain, imm_cdr, imm_split, imm_tcrnum           # noqa: E402


def atom(serial, name, resn, chain, resseq, icode=" ", x=0.0, rec="ATOM  ", el="C"):
    return (f"{rec}{serial:5d} {name:<4s} {resn:3s} {chain}{resseq:4d}{icode}   "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{20.0:6.2f}          {el:>2s}\n")


def chain_lines(chain, n, start=1, x0=0.0):
    """n alanines (CA only) 3.8 A apart along x."""
    return [atom(i + 1, "CA", "ALA", chain, start + i, x=x0 + 3.8 * i) for i in range(n)]


# ------------------------------------------------------------------ registry
def test_registry_is_consistent():
    assert set(_registry.ROLE_OF_KIND.values()) <= set(_registry.ROLES)
    assert set(_registry.LAYOUT) == set(_registry.ROLES)
    ids = list(_registry.LAYOUT.values())
    assert len(ids) == len(set(ids)) and all(len(c) == 1 for c in ids)
    assert len({c.upper() for c in ids}) == len(ids)       # no clash for case-insensitive viewers
    for kind in list(_registry.ANARCI_KIND.values()) + list(_registry.MHC_KIND.values()):
        assert kind in _registry.KINDS
    for _, roles in _registry.PARTS.values():
        assert set(roles) <= set(_registry.ROLES)
    for kind in _registry.CDR_SIDE:
        assert _registry.family(kind) in ("tcr", "antibody")


# ------------------------------------------------------------------ _io
def test_check_input_options(tmp_path):
    f = tmp_path / "x.pdb"
    f.write_text("ATOM\n")
    opt, fh = _io.check_input(["-A,B", str(f)], "doc")
    assert opt == "A,B" and fh.name == str(f)
    fh.close()
    opt, fh = _io.check_input([str(f)], "doc")
    assert opt == "" and fh.name == str(f)
    fh.close()
    with pytest.raises(SystemExit):
        _io.check_input(["A", str(f)], "doc")                # option without '-'
    with pytest.raises(SystemExit):
        _io.check_input([str(tmp_path / "missing.pdb")], "doc")
    cif = tmp_path / "x.cif"
    cif.write_text("data_x\n")
    opt, fh = _io.check_input([str(cif)], "doc")            # mmCIF is accepted
    fh.close()


def test_chain_list():
    assert _io.chain_list("") is None
    assert _io.chain_list("A,B") == {"A", "B"}
    assert _io.chain_list("AA,B1") == {"AA", "B1"}            # mmCIF chain IDs


def test_apply_numbering_changes_only_number_columns():
    lines = chain_lines("A", 3) + ["TER       4      ALA A   3\n"] + chain_lines("B", 2)
    new = {("A", 1, "", False): (10, ""), ("A", 2, "", False): (10, "A"),
           ("A", 3, "", False): (11, "")}
    out = _io.apply_numbering(lines, new)
    assert [l[22:27] for l in out[:3]] == ["  10 ", "  10A", "  11 "]
    assert out[3][22:27] == "  11 "                            # TER follows its chain
    assert out[4:] == lines[4:]                                # chain B untouched
    for a, b in zip(lines[:3], out[:3]):
        assert a[:22] == b[:22] and a[27:] == b[27:]


def test_apply_numbering_moves_colliding_water():
    lines = chain_lines("A", 2) + [atom(9, "O", "HOH", "A", 100, rec="HETATM", el="O")]
    new = {("A", 1, "", False): (100, ""), ("A", 2, "", False): (101, "")}
    out = _io.apply_numbering(lines, new)
    assert out[2][22:26] == " 102"                             # moved after the chain


def test_header_records_follow_and_ss_ranges_split():
    # residues 1..6 become 1, 2, 2F, 2E (descending insertion codes, as in IMGT), 3, 4
    lines = ["SHEET    1   S 1 ALA A   1  ALA A   6  0\n",
             "SSBOND   1 CYS A    2    CYS A    5                          1555   1555  2.03\n"]
    lines += chain_lines("A", 6)
    new = {("A", i + 1, "", False): v for i, v in
           enumerate([(1, ""), (2, ""), (2, "F"), (2, "E"), (3, ""), (4, "")])}
    out = _io.apply_numbering(lines, new)
    sheets = [l for l in out if l.startswith("SHEET")]
    rng = [((r[0][1], r[0][2]), (r[1][1], r[1][2])) for r in
           ([x[:3] for x in _io._refs(s)] for s in sheets)]
    # every range ascends and no range captures a residue of another range
    keys = [(1, ""), (2, ""), (2, "F"), (2, "E"), (3, ""), (4, "")]
    covered = []
    for lo, hi in rng:
        assert lo <= hi
        covered += [k for k in keys if lo <= k <= hi]
    assert sorted(covered) == sorted(keys)                     # each residue exactly once
    ss = [l for l in out if l.startswith("SSBOND")][0]
    assert [r[1:3] for r in _io._refs(ss)] == [(2, ""), (3, "")]


def test_rename_and_select_chains():
    lines = ["SSBOND   1 CYS A    1    CYS B    1                          1555   1555  2.03\n"]
    lines += chain_lines("A", 2) + chain_lines("B", 2, x0=50)
    ren = _io.rename_chains(lines, {"A": "M", "B": "2"})
    assert {l[21] for l in ren if l.startswith("ATOM")} == {"M", "2"}
    assert [r[0] for r in _io._refs(ren[0])] == ["M", "2"]
    sel = _io.select_chains(lines, {"A"})
    assert {l[21] for l in sel if l.startswith("ATOM")} == {"A"}
    assert not any(l.startswith("SSBOND") for l in sel)       # refers to a removed chain


# ------------------------------------------------------------------ _mhc
def test_fill_symmetric_loop():
    slots = [(84, "A"), (84, "B"), (84, "C"), (85, "C"), (85, "B"), (85, "A")]
    assert _mhc._fill(slots, 3, (84, "")) == [(84, "A"), (84, "B"), (85, "A")]
    assert _mhc._fill(slots, 2, (84, "")) == [(84, "A"), (85, "A")]


def test_fill_more_residues_than_slots():
    assert _mhc._fill([(61, "A")], 3, (61, "")) == [(61, "A"), (61, "B"), (61, "C")]


def test_mhc_numbering_hla_a2_sequence():
    hla_a2 = ("GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDL"
              "GTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQL"
              "RAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRP"
              "AGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWE")
    labels, info = _mhc.number_sequence(hla_a2)
    assert info["chain_type"] == "MHC_I_ALPHA"
    assert labels[0] == (1, "") and labels[100] == (1011, "")   # alpha2 Cys (HLA-A2 C101)
    assert hla_a2[100] == "C"
    assert len(set(labels)) == len(labels)


def test_tcr_constant_domain_is_not_mhc():
    trbc = ("EDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVSTDPQPLKEQPALNDSRYALSSRL"
            "RVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD")
    with pytest.raises(ValueError):
        _mhc.number_sequence(trbc)


# ------------------------------------------------------------------ _tcr
def test_continue_after():
    labels = [None, (1, ""), (2, ""), (128, ""), None, None]
    assert _tcr.continue_after(labels) == [(0, ""), (1, ""), (2, ""), (128, ""), (129, ""), (130, "")]


@pytest.mark.skipif(not _tcr.available(), reason="ANARCI not installed")
def test_tcr_full_labels_with_constant():
    seq = ("DAGVIQSPRHEVTEMGQEVTLRCKPISGHNSLFWYRQTMMRGLELLIYFNNNVPIDDSGMPEDRFSAKMPNASFSTLKIQ"
           "PSEPRDSAVYFCASSYVGNTGELFFGEGSRLTVLEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSW"
           "WVNGKEVHSGVSTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEA"
           "WGRAD")
    labels, info = _tcr.full_labels(seq, constant=True)
    assert info["domains"][0]["chain_type"] == "B" and info["constant"][0] == "C_TRBC"
    assert len(set(labels)) == len(labels)
    num = {l: seq[i] for i, l in enumerate(labels)}
    assert num[(23, "")] == "C" and num[(104, "")] == "C"           # V-domain anchors
    assert num[(1023, "")] == "C" and num[(1104, "")] == "C"         # C-domain anchors


# ------------------------------------------------------------------ _annot
def test_identity():
    assert _annot._identity("ACDEFG", "ACDEFG") == 1.0
    assert _annot._identity("ACDEFG", "ACDQFG") == pytest.approx(5 / 6)


# ------------------------------------------------------------------ _complex
class FakeModel(list):
    """Stand-in for a Bio.PDB model: a list of chains with .id."""


def complex_ctx(kinds, contacts):
    from Bio.PDB.Chain import Chain
    model = FakeModel(Chain(c) for c in kinds)
    ident = {c: {"kind": k} for c, k in kinds.items()}
    ctx = _complex.Context.__new__(_complex.Context)
    ctx.ident = ident
    ctx.role = {c: _registry.ROLE_OF_KIND.get(k) for c, k in kinds.items()}
    ctx.ids = [c for c in kinds if ctx.role[c]]
    ctx.all = list(kinds)
    ctx.C = {tuple(sorted(k)): v for k, v in contacts.items()}
    ctx.free = set(ctx.all)
    return model, ctx


def test_tcr_pmhc_assembler_two_copies():
    kinds = {"A": "MHC_I", "B": "B2M", "C": "peptide", "D": "TCR_alpha", "E": "TCR_beta",
             "F": "MHC_I", "G": "B2M", "H": "peptide", "I": "TCR_alpha", "J": "TCR_beta"}
    con = {("A", "B"): 50, ("A", "C"): 30, ("A", "D"): 20, ("C", "D"): 5, ("D", "E"): 60,
           ("F", "G"): 50, ("F", "H"): 30, ("F", "I"): 20, ("I", "J"): 60,
           ("A", "G"): 2}                                           # crystal contact, too few
    _, ctx = complex_ctx(kinds, con)
    cxs = _complex.tcr_pmhc(ctx)
    assert cxs == [{"mhc": "A", "b2m": "B", "pep": "C", "tcra": "D", "tcrb": "E"},
                   {"mhc": "F", "b2m": "G", "pep": "H", "tcra": "I", "tcrb": "J"}]
    assert not ctx.free & set("ABCDEFGHIJ")


def test_tcr_only_assembler():
    _, ctx = complex_ctx({"D": "TCR_alpha", "E": "TCR_beta", "X": "TCR_beta"},
                         {("D", "E"): 40})
    assert _complex.tcr_only(ctx) == [{"tcra": "D", "tcrb": "E"}, {"tcrb": "X"}]


# ------------------------------------------------------------------ tool options
def test_rechain_layout_option():
    assert imm_rechain.layout("")["b2m"] == "2"
    lay = imm_rechain.layout("tcra:D,tcrb:E")
    assert lay["tcra"] == "D" and lay["tcrb"] == "E" and lay["mhc"] == "M"
    lay = imm_rechain.layout("D,E,A,X,B,C")
    assert [lay[r] for r in _registry.ROLES[:6]] == list("DEAXBC")
    with pytest.raises(SystemExit):
        imm_rechain.layout("tcra:DD")


def test_cdr_selection_option():
    assert imm_cdr.selection("cdr3") == {("a", 3), ("b", 3)}
    assert imm_cdr.selection("a1,b2") == {("a", 1), ("b", 2)}
    with pytest.raises(SystemExit):
        imm_cdr.selection("x9")


def test_split_parts_and_tcrnum_option():
    assert imm_split.parts_of("") == list(_registry.PARTS)
    assert imm_split.parts_of("tcr,pep") == ["tcr", "pep"]
    with pytest.raises(SystemExit):
        imm_split.parts_of("fab")
    assert imm_tcrnum.parse("const") == (True, None)
    assert imm_tcrnum.parse("const:D,E") == (True, {"D", "E"})
    assert imm_tcrnum.parse("D") == (False, {"D"})


# ------------------------------------------------------------------ _doc (PDB / mmCIF)
DATA = Path(__file__).parent / "data"


def test_read_detects_format():
    from immtools import _doc
    assert _doc.read(open(DATA / "3qdg_rcsb.cif")).fmt == "cif"
    assert _doc.read(open(DATA / "3qdg_rcsb.pdb")).fmt == "pdb"


def test_cif_and_pdb_give_the_same_numbering():
    from immtools import _doc, _number
    res = []
    for f in ("3qdg_rcsb.pdb", "3qdg_rcsb.cif"):
        d = _doc.read(open(DATA / f))
        new = _number.mhc_ids(d.model(), quiet=True)
        model = d.renumbered(new).model()
        res.append([(c.id, [r.id[1:] for r in c]) for c in model])
    assert res[0] == res[1]


def test_pdb_rejects_long_chain_ids():
    from immtools import _doc
    d = _doc.read(open(DATA / "3qdg_rcsb.pdb"))
    with pytest.raises(SystemExit):
        d.renamed({"A": "MHC"})
