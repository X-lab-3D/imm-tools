"""Write examples/try_imm_tools.ipynb (the Colab notebook linked from the README).

    python dev/make_notebook.py [test.ipynb]

With a file name, also writes a copy without the install and upload cells, to run locally:
    jupyter nbconvert --to notebook --execute test.ipynb
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import nbformat as nbf
md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
cells = [
md("""# Try imm-tools

Consistent IMGT numbering and chain IDs for TCR–pMHC structures. Click a cell and press
the ▶ button (or **Runtime → Run all**). Everything runs in this free Colab session; nothing
is installed on your computer.

**Lines starting with `!` are terminal commands.** On your own computer, after installing
imm-tools, type the same commands in a terminal **without the `!`**, e.g.
`imm_annotate 3qdg.pdb`.

Documentation: <https://github.com/X-lab-3D/imm-tools>"""),
md("## 1. Install (about a minute)"),
code("""!pip install -q git+https://github.com/X-lab-3D/imm-tools.git py3Dmol
!imm_setup"""),
md("""## 2. Get a structure

3QDG: the DMF5 TCR bound to HLA-A2 presenting a MART-1 peptide. (To try another entry,
replace `3QDG` / `3qdg` in the cells below.)"""),
code("""!curl -s -o 3qdg.pdb https://files.rcsb.org/download/3QDG.pdb"""),
md("""## 3. What is in the file?

The output is a tab-separated table: one row per chain."""),
code("""!imm_annotate 3qdg.pdb"""),
md("""## 4. One complex, standard chain IDs, IMGT numbering

Each tool does one step and reports what it did. `imm_complex -1` keeps the first complex
(some files hold several copies)."""),
code("""!imm_complex -1 3qdg.pdb > step1.pdb
!imm_rechain step1.pdb > step2.pdb
!imm_tcrnum step2.pdb > step3.pdb
!imm_mhcnum step3.pdb > step4.pdb
!imm_pepnum step4.pdb > 3qdg_imgt.pdb"""),
md("The same in one line: `|` passes each tool's output to the next, so no intermediate files are needed."),
code("""!imm_complex -1 3qdg.pdb | imm_rechain | imm_tcrnum | imm_mhcnum | imm_pepnum > 3qdg_imgt.pdb"""),
md("""Now the structure has TCR α = chain **A**, TCR β = **B**, MHC = **M**, β2m = **2**,
peptide = **P**, and IMGT residue numbers – like every other structure you treat this way.
Check it:"""),
code("""!imm_check 3qdg_imgt.pdb"""),
md("## 5. CDR loops"),
code("""!imm_cdr 3qdg_imgt.pdb"""),
md("""## 6. Closest common allele (human TCRs)

`imm_annotate` also reports, for the V and J gene, the closest allele among those common in
the population (frequency ≥ 5%; default source: VDJ.online, TCR repertoires of mostly European
donors) and how the chain differs from it. For an engineered TCR the differences are the
mutations. 5BRZ is an affinity-enhanced TCR:"""),
code("""!curl -s -o 5brz.pdb https://files.rcsb.org/download/5BRZ.pdb
!imm_annotate -D,E 5brz.pdb | cut -f2,7,9,11,12"""),
md("""Choose another frequency source with `-freq:<source>` (and a population with
`-freq:<source>/<population>`); the README lists the available ones."""),
md("""## 7. Look at it

(This cell is Python, for the 3D view in the notebook.) Because chains and numbers are
standard, the same lines colour **any** TCR–pMHC structure: CDR1 yellow, CDR2 green, CDR3
orange, groove helices blue, peptide pink."""),
code("""import py3Dmol

view = py3Dmol.view(width=700, height=500)
view.addModel(open("3qdg_imgt.pdb").read(), "pdb")
view.setStyle({"cartoon": {"color": "#e1e0d9"}})
view.setStyle({"chain": "M", "resi": ["50-92", "1050-1092"]}, {"cartoon": {"color": "#2a78d6"}})
view.setStyle({"chain": "N", "resi": "50-92"}, {"cartoon": {"color": "#2a78d6"}})    # class II
view.setStyle({"chain": ["A", "B"], "resi": "27-38"}, {"cartoon": {"color": "#eda100"}})
view.setStyle({"chain": ["A", "B"], "resi": "56-65"}, {"cartoon": {"color": "#1baf7a"}})
view.setStyle({"chain": ["A", "B"], "resi": "105-117"}, {"cartoon": {"color": "#eb6834"}})
view.setStyle({"chain": "P"}, {"stick": {"color": "#e87ba4"}})
view.zoomTo()
view.show()"""),
md("""## 8. Your own files

(Python, for uploading in Colab.) Upload PDB or mmCIF files – they stay in this session –
and download them standardised:"""),
code("""from google.colab import files
for name in files.upload():
    stem, ext = name.rsplit(".", 1)
    !imm_complex -1 "{name}" | imm_rechain | imm_tcrnum | imm_mhcnum | imm_pepnum > "{stem}_imgt.{ext}"
    files.download(f"{stem}_imgt.{ext}")"""),
md("""## More

* Help for any tool: `!imm_rechain -h` (in a terminal: `imm_rechain -h`).
* More examples (whole folders, docking RMSD, contact maps, PyMOL) and the Python interface:
  [README](https://github.com/X-lab-3D/imm-tools) and
  [examples/](https://github.com/X-lab-3D/imm-tools/tree/main/examples)."""),
]
nb = nbf.v4.new_notebook(cells=cells, metadata={"kernelspec": {"name": "python3", "display_name": "Python 3"},
                                                "colab": {"provenance": []}})
nbf.write(nb, ROOT / "examples" / "try_imm_tools.ipynb")
t = nbf.v4.new_notebook(cells=[c for c in cells if "pip install" not in c.source and "google.colab" not in c.source],
                        metadata=nb.metadata)
if len(sys.argv) > 1:      # test copy without the Colab-only cells, to run locally
    nbf.write(t, sys.argv[1])
