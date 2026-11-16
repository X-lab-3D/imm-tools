# dev

Scripts that built some of the shipped data (not installed; they need the source data):

* `build_tcr_constant.py` - TCR constant-domain HMMs (`immtools/models/tcr_constant`) from IMGT/DomainDisplay
* `validate_tcr_constant.py` - compare with IMGT/3Dstructure-DB numbering
* `build_mhc_refdb.py` - MHC reference sequences for imm_annotate (`immtools/models/mhc_refs.fasta.gz`)
* `make_notebook.py` - writes the Colab notebook `examples/try_imm_tools.ipynb`
* `alleles/` - allele-frequency sources for the closest common TCR allele
  (`immtools/models/tcr_alleles/<source>.tsv`): `germline.py` (shared numbering and table
  writing), one `build_<source>.py` per source. `build_vdjonline.py` reads
  `data/vdjonline/` (CSV downloads of https://vdj.online/library, 2026-11-16).

The MHC HMMs (`immtools/models/mhc_broad`, `mhc_mammal`) were built with a separate pipeline
that is not part of this repository.
