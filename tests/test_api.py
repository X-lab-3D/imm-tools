"""Tests of the public Python interface (import immtools)."""
from pathlib import Path

import pytest

import immtools

DATA = Path(__file__).parent / "data"
HLA_A2 = ("GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYN"
          "QSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWL"
          "RRYLENGK")
TRB = ("DAGVIQSPRHEVTEMGQEVTLRCKPISGHNSLFWYRQTMMRGLELLIYFNNNVPIDDSGMPEDRFSAKMPNASFSTLKIQPSEPR"
       "DSAVYFCASSYVGNTGELFFGEGSRLTVLEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSG"
       "VSTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD")


def test_number_mhc():
    n = immtools.number_mhc(HLA_A2)
    assert n.chain_type == "MHC_I"
    assert n.labels[:3] == ["1", "2", "3"] and dict(zip(n.labels, n.sequence))["1011"] == "C"
    assert [d[0] for d in n.domains] == ["G_ALPHA1", "G_ALPHA2"]
    with pytest.raises(ValueError):
        immtools.number_mhc("ACDEFGHIKL" * 10)


def test_number_tcr():
    n = immtools.number_tcr(TRB, constant=True)
    lab = dict(zip(n.labels, n.sequence))
    assert n.chain_type == "TCR_beta" and n.genes["v_gene"].startswith("TRBV")
    assert lab["23"] == lab["104"] == lab["1023"] == lab["1104"] == "C"
    with pytest.raises(ValueError):
        immtools.number_tcr(HLA_A2)


def test_structure_functions(tmp_path):
    pdb = DATA / "3qdg_rcsb.pdb"
    rows = immtools.annotate(pdb)
    assert [r["type"] for r in rows] == ["MHC_I", "B2M", "peptide", "TCR_alpha", "TCR_beta"]
    assert rows[0]["gene_allele"].startswith("HLA-A*02:01")
    assert immtools.complexes(DATA / "1d9k_two_copies.pdb")[1] == \
        {"mhc": "G", "mhcb": "H", "pep": "Q", "tcra": "E", "tcrb": "F"}
    cdr3 = {c["chain"]: c["sequence"] for c in immtools.cdrs(pdb) if c["cdr"] == "CDR3"}
    assert set(cdr3) == {"D", "E"}
    out = tmp_path / "x.cif"
    text = immtools.renumber(DATA / "3qdg_rcsb.cif", out=out, constant=True)
    assert text.startswith("data_") and out.read_text() == text
    assert {r["status"] for r in immtools.check(text) if r["check"] == "numbering"} == {"OK"}
    std = immtools.rechain(immtools.select_complex(DATA / "1d9k_two_copies.pdb", 2))
    assert [c["chain"] for c in immtools.annotate(std)] == ["A", "B", "M", "N", "P"]
    with pytest.raises(ValueError):
        immtools.rechain(DATA / "1d9k_two_copies.pdb")       # two complexes
    with open(pdb) as fh:                                    # open files work too
        assert len(immtools.annotate(fh)) == 5
    model = immtools.structure(text)                         # Biopython model of the result
    assert model["D"][104].get_resname() == "CYS"
    assert model["A"][2023].get_resname() == "CYS"


def test_examples(tmp_path):
    """The example scripts run on the test structures."""
    import subprocess
    import sys
    ex = Path(__file__).parents[1] / "examples"
    r = subprocess.run([sys.executable, ex / "cdr3_table.py", DATA / "1mi5_tcr_pmhc.pdb"],
                       capture_output=True, text=True, check=True)
    rows = [l.split("\t") for l in r.stdout.splitlines()]
    assert rows[1][5] == "ILPLAGGTSYGKLT" and "112A" in rows[1][7].split(",")
    imgt = tmp_path / "3qdg_imgt.pdb"
    immtools.renumber(DATA / "3qdg_rcsb.pdb", out=imgt)
    r = subprocess.run([sys.executable, ex / "pymol_regions.py", "-footprint", imgt],
                       capture_output=True, text=True, check=True)
    assert "select cdr3b, 3qdg_imgt and chain E and resi 105+106" in r.stdout
    assert "python end" in r.stdout
