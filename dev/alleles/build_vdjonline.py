#!/usr/bin/env python3
"""Allele-frequency source 'vdjonline': human TCR alleles from the VDJ.online gene library.

Input: CSV downloads of https://vdj.online/library (MiLaboratories), one per chain and segment,
in dev/data/vdjonline/ (downloaded 2026-11-16). Only human alleles with a frequency are used;
frequencies come from immune repertoire sequencing of 134 donors of mostly European descent
(Mikelov et al. 2024). VDJ.online places no restrictions on the use or distribution of its
molecular data.

    python dev/alleles/build_vdjonline.py    ->   immtools/models/tcr_alleles/vdjonline.tsv
"""
import csv
from pathlib import Path

from germline import ROOT, number, translations, write_source

SRC = ROOT / "dev" / "data" / "vdjonline"
META = dict(
    name="vdjonline",
    title="VDJ.online gene library (MiLaboratories)",
    species="human",
    populations="all = TCR repertoires of 134 donors, predominantly of European descent "
                "(about 210-260 haplotypes per gene)",
    citation="A. Mikelov et al., Ultrasensitive allele inference from immune repertoire "
             "sequencing data with MiXCR, Genome Res. 34, 2293-2303 (2024), "
             "https://pmc.ncbi.nlm.nih.gov/articles/PMC11694755/; MiXCR: D. A. Bolotin et al., "
             "Nat. Methods 12, 380-381 (2015)",
    url="https://vdj.online/library",
    licence="no restrictions on use or distribution of the molecular data (VDJ.online terms); "
            "cite MiXCR",
)


def rows():
    left_out = []
    for f in sorted(SRC.glob("*.csv")):
        for r in csv.DictReader(open(f)):
            if r["species"] != "hs" or float(r["freq"]) <= 0 or r["geneType"] not in ("V", "J"):
                continue
            seg = r["geneType"]
            protein = r["v-region aa"] if seg == "V" else r["aa"]
            nts = r["v-region nt"] if seg == "V" else r["nt"]
            imgt = next((x for x in (number(p, r["geneName"], seg)
                                     for p in translations(protein, nts)) if x), None)
            if imgt is None:
                left_out.append(r["name"])
                continue
            yield dict(gene=r["geneName"].replace("/", ""), allele=r["name"].replace("/", ""),
                       segment=seg, population="all", freq=f"{float(r['freq']):.4f}",
                       n=r["supportingHaplotypes"], imgt=imgt)
    print("left out (stop codons in every frame, i.e. not functional, or no ANARCI germline "
          f"of the gene): {', '.join(left_out)}")


if __name__ == "__main__":
    path, n = write_source(META, list(rows()), __file__)
    print(f"{path}: {n} alleles")
