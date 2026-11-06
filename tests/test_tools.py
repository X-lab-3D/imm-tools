"""Run with:  pytest tests/"""
import os
import subprocess
from pathlib import Path

DATA = Path(__file__).parent / "data"
PDB1 = DATA / "1mi5_tcr_pmhc.pdb"       # HLA-B*08:01 + b2m + peptide + TCR (LC13)
PDB2 = DATA / "2pxy_tcr_pmhc2.pdb"      # mouse I-A(u) + TCR
TWO = DATA / "1d9k_two_copies.pdb"      # two TCR-pMHC-II complexes in the asymmetric unit


def run(*args, inp=None, cwd=None):
    p = subprocess.run([str(a) for a in args], input=inp, capture_output=True, text=True, cwd=cwd)
    assert p.returncode == 0, p.stderr
    return p.stdout, p.stderr


def ca(text, chain):
    return [l[17:20] + l[22:27].strip() for l in text.splitlines()
            if l.startswith("ATOM") and l[12:16] == " CA " and l[21] == chain]


def tsv(text):
    rows = [l.split("\t") for l in text.strip().splitlines()]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def test_seq_types():
    out, _ = run("imm_seq", PDB1)
    assert [h.split("|")[1] for h in out.splitlines() if h.startswith(">")] == \
        ["MHC_I", "B2M", "peptide", "TCR_alpha", "TCR_beta"]
    out, _ = run("imm_seq", "-mhc", PDB2)
    assert [h.split("|")[1] for h in out.splitlines() if h.startswith(">")] == \
        ["MHC_II_alpha", "MHC_II_beta"]


def test_annotate():
    out, _ = run("imm_annotate", PDB1)
    rows = {r["chain"]: r for r in tsv(out)}
    assert rows["A"]["gene_allele"].startswith("HLA-B*08:01") and rows["A"]["identity"] == "100.0"
    assert rows["D"]["gene_allele"].startswith("TRAV26-2") and rows["E"]["gene_allele"].startswith("TRBV7-8")
    assert rows["A"]["domains"] == "G_ALPHA1:1-90;G_ALPHA2:91-182;C_ALPHA3:183-277"
    assert {r["complex"] for r in rows.values()} == {"1"} and "role" not in rows["A"]


def test_mhcnum_class_I():
    out, err = run("imm_mhcnum", PDB1)
    a = ca(out, "A")
    assert a[:3] == ["GLY1", "SER2", "HIS3"]
    assert "ALA1092" in a and "ASP2001A" in a           # IMGT alpha2/alpha3 junction
    assert ca(out, "B")[:4] == ["ILE1C", "GLN1B", "ARG1A", "THR1"]
    assert ca(out, "D") == ca(PDB1.read_text(), "D")     # TCR untouched
    assert "left unchanged" in err


def test_mhcnum_class_II():
    out, _ = run("imm_mhcnum", PDB2)
    b = ca(out, "B")
    assert "GLN61A" in b and "GLU1001A" in b and "GLN1001" in b


def test_only_number_columns_change():
    src = PDB1.read_text()
    out, _ = run("imm_mhcnum", PDB1)
    out, _ = run("imm_tcrnum", "-const", inp=out)
    out, _ = run("imm_pepnum", inp=out)
    a, b = src.splitlines(), out.splitlines()
    assert len(a) == len(b)
    for x, y in zip(a, b):
        if x.startswith(("ATOM", "HETATM")):
            assert x[:22] == y[:22] and x[27:] == y[27:]


def test_stdin_and_chain_option():
    out1, _ = run("imm_mhcnum", "-A", PDB1)
    out2, _ = run("imm_mhcnum", "-A", inp=PDB1.read_text())
    assert out1 == out2
    assert ca(out1, "B") == ca(PDB1.read_text(), "B")


def test_tcrnum_v_and_constant():
    out, _ = run("imm_tcrnum", PDB1)
    for ch in "DE":
        res = ca(out, ch)
        assert "CYS23" in res and "CYS104" in res and "TRP41" in res
    out, _ = run("imm_tcrnum", "-const", PDB1)
    e = ca(out, "E")
    assert "CYS1023" in e and "CYS1104" in e and "TRP1041" in e      # TRBC C-LIKE anchors
    assert ca(out, "A") == ca(PDB1.read_text(), "A")


def test_pepnum():
    src = PDB1.read_text().replace(" C   1 ", " C 101 ")
    out, _ = run("imm_pepnum", inp=src)
    assert [x[3:] for x in ca(out, "C")] == [str(i) for i in range(1, 10)]


def test_check_pipeline():
    out, _ = run("imm_mhcnum", PDB1)
    out, _ = run("imm_tcrnum", inp=out)
    out, _ = run("imm_check", inp=out)
    rows = tsv(out)
    num = {r["chain"]: r["status"] for r in rows if r["check"] == "numbering"}
    anc = {r["chain"]: r["status"] for r in rows if r["check"] == "anchors"}
    assert set(num.values()) == {"OK"} and set(anc.values()) == {"OK"}
    raw = tsv(run("imm_check", PDB1)[0])
    assert any(r["chain"] == "A" and r["check"] == "numbering" and r["status"] == "FAIL" for r in raw)


def test_cdr():
    out, _ = run("imm_cdr", PDB1)
    cdr3 = {r["chain"]: r["sequence"] for r in tsv(out) if r["cdr"] == "CDR3"}
    assert cdr3 == {"D": "ILPLAGGTSYGKLT", "E": "ASSLGQAYEQY"}
    out, _ = run("imm_cdr", "-b3", PDB1)
    assert len(ca(out, "E")) == 11 and not ca(out, "D")


def test_complex_rechain_split(tmp_path):
    out, err = run("imm_complex", "-2", TWO)
    assert "complex 2: tcra=E tcrb=F mhc=G mhcb=H pep=Q" in err
    assert {l[21] for l in out.splitlines() if l.startswith("ATOM")} == set("EFGHQ")
    std, _ = run("imm_rechain", inp=out)
    assert [c for c in dict.fromkeys(l[21] for l in std.splitlines() if l.startswith("ATOM"))] == \
        list("ABMNP")
    out, _ = run("imm_complex", "-1", DATA / "3qdg_rcsb.pdb")
    std, _ = run("imm_rechain", inp=out)
    assert "".join(dict.fromkeys(l[21] for l in std.splitlines() if l.startswith("ATOM"))) == "M2PAB"
    p = subprocess.run(["imm_rechain", str(TWO)], capture_output=True, text=True)
    assert p.returncode == 1 and "2 complexes" in p.stderr
    run("imm_split", "-tcr,pmhc", PDB1.resolve(), cwd=tmp_path)
    assert sorted(os.listdir(tmp_path)) == ["1mi5_tcr_pmhc_TCR.pdb", "1mi5_tcr_pmhc_pMHC.pdb"]


def test_water_collision():
    src = PDB1.read_text().splitlines(keepends=True)
    water = "HETATM 9999  O   HOH A1001      10.000  10.000  10.000  1.00 20.00           O  \n"
    i = max(k for k, l in enumerate(src) if l.startswith("ATOM") and l[21] == "A")
    src.insert(i + 1, water)
    out, _ = run("imm_mhcnum", "-A", inp="".join(src))
    w = [l for l in out.splitlines() if l.startswith("HETATM") and l[17:20] == "HOH"][0]
    nums = {l[22:27] for l in out.splitlines() if l.startswith("ATOM") and l[21] == "A"}
    assert w[22:27] not in nums


def test_header_records_follow_numbering():
    """HELIX / SHEET / SSBOND point at residues that exist after renumbering (PyMOL reads them)."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from immtools import _io
    out, _ = run("imm_tcrnum", "-const", DATA / "3qdg_rcsb.pdb")
    out, _ = run("imm_mhcnum", inp=out)
    out, _ = run("imm_rechain", inp=out)
    lines = out.splitlines(keepends=True)
    atoms = {(l[21], int(l[22:26]), l[26].strip()) for l in lines if l.startswith(("ATOM", "HETATM"))}
    refs = [r[:3] for l in lines if l[:6] in ("HELIX ", "SHEET ", "SSBOND") for r in _io._refs(l)]
    assert refs and all(r in atoms for r in refs)
    assert any(r[0] == "B" and r[1] > 1000 for r in refs)       # TCR beta constant domain


def _ss_by_residue(text):
    """Secondary structure per residue (file order) the way PyMOL reads HELIX / SHEET:
    every residue whose (number, insertion code) lies between start and end."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from immtools import _io
    lines = text.splitlines(keepends=True)
    res = list(dict.fromkeys((l[21], int(l[22:26]), l[26].strip()) for l in lines
                             if l.startswith("ATOM") and l[12:16] == " CA "))
    ss = {r: "L" for r in res}
    for l in lines:
        if l[:6] in ("HELIX ", "SHEET "):
            refs = _io._refs(l)
            (c, n0, i0, _), (_, n1, i1, _) = refs[0], refs[1]
            for r in res:
                if r[0] == c and (n0, i0) <= (r[1], r[2]) <= (n1, i1):
                    ss[r] = l[0]
    return [ss[r] for r in res]


def test_secondary_structure_unchanged():
    """After renumbering (descending IMGT insertion codes!) viewers see the same SS."""
    src = (DATA / "3qdg_rcsb.pdb").read_text()
    out, _ = run("imm_tcrnum", "-const", DATA / "3qdg_rcsb.pdb")
    out, _ = run("imm_mhcnum", inp=out)
    assert _ss_by_residue(out) == _ss_by_residue(src)



def test_mmcif_in_mmcif_out(tmp_path):
    """mmCIF input gives mmCIF output with the same numbering as the PDB route."""
    import gemmi
    cif = DATA / "3qdg_rcsb.cif"
    out, _ = run("imm_tcrnum", "-const", cif)
    out, _ = run("imm_mhcnum", inp=out)
    assert out.startswith("data_")
    pdb_out, _ = run("imm_tcrnum", "-const", DATA / "3qdg_rcsb.pdb")
    pdb_out, _ = run("imm_mhcnum", inp=pdb_out)
    (tmp_path / "x.cif").write_text(out)
    a = gemmi.read_structure(str(tmp_path / "x.cif"))
    b = gemmi.read_pdb_string(pdb_out)
    def ids(st):
        out = {}
        for ch in st[0]:
            out.setdefault(ch.name, []).extend(str(r.seqid) for r in ch if r.het_flag == "A")
        return out
    assert ids(a) == ids(b)
    rows = tsv(run("imm_annotate", cif)[0])
    assert [r["type"] for r in rows] == ["MHC_I", "B2M", "peptide", "TCR_alpha", "TCR_beta"]
    run("imm_split", "-tcr", cif.resolve(), cwd=tmp_path)
    assert (tmp_path / "3qdg_rcsb_TCR.cif").exists()


def test_mmcif_long_chain_ids(tmp_path):
    """Chain IDs longer than one character (only possible in mmCIF)."""
    import gemmi
    doc = gemmi.cif.read(str(DATA / "3qdg_rcsb.cif"))
    b = doc.sole_block()
    at = b.get_mmcif_category("_atom_site.")
    at["auth_asym_id"] = [c + "x" for c in at["auth_asym_id"]]
    b.set_mmcif_category("_atom_site.", at)
    f = tmp_path / "long.cif"
    doc.write_file(str(f))
    out, err = run("imm_rechain", f)
    assert "tcra: Dx -> A" in err and "mhc: Ax -> M" in err
    out, _ = run("imm_mhcnum", "-Ax", f)
    st = gemmi.read_structure_string(out) if hasattr(gemmi, "read_structure_string") else None
    (tmp_path / "n.cif").write_text(out)
    st = gemmi.read_structure(str(tmp_path / "n.cif"))
    nums = [str(r.seqid) for r in st[0]["Ax"] if r.het_flag == "A"]
    assert "1092" in nums and "2001A" in nums
