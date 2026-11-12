# Benchmarks

Accuracy of the numbering, measured as the fraction of residues that get the expected IMGT
number. Model version: see `immtools/models/MODELS_VERSION.json`. "Broad" and "mammal" are
the two MHC model sets; "combined" is what the tools use (best-scoring set per domain).

## MHC: human and mouse structures vs IMGT/3Dstructure-DB

All MHC chains of 78 IMGT/3Dstructure-DB entries (human and mouse class I and II), compared
with IMGT's own residue numbering. **These structures were part of the training data**, so
this shows that the models reproduce IMGT, not how well they generalise.

| chain | chains | residues | broad | combined |
|---|---|---|---|---|
| class I heavy chain | 78 | 21,162 | 100.00 % | 100.00 % |
| β2m | 79 | 7,806 | 100.00 % | 100.00 % |
| class II α | 40 | 7,076 | 99.96 % | 99.96 % |
| class II β | 23 | 4,174 | 98.99 % | 99.74 % |

Remaining differences: mouse I-A β position 61A vs IMGT 62 (the structure supports 61A), one
class II α junction residue, and residues at the end of the β1 helix of H2-ABu.

## MHC: species left out

For each species with experimental structures, all of its sequences and structures were
removed from training, the models rebuilt, and that species' structures numbered. Reference:
IMGT numbering after the structure check (see `differences_from_imgt.tsv`); "convention"
also applies the loop convention used here where IMGT numbers the same loop length
differently in different species.

| species left out | chain | mammal | broad | combined | combined (convention) |
|---|---|---|---|---|---|
| *Bos taurus* | class I | 100.00 % | 100.00 % | 100.00 % | 100.00 % |
| | β2m | 99.62 % | 100.00 % | 99.62 % | 99.62 % |
| *Ctenopharyngodon idella* (grass carp) | class I | 93.61 % | 91.03 % | 93.61 % | 93.61 % |
| | β2m | 100.00 % | 98.96 % | 100.00 % | 100.00 % |
| *Equus caballus* | class I | 100.00 % | 100.00 % | 100.00 % | 100.00 % |
| | β2m | 100.00 % | 100.00 % | 100.00 % | 100.00 % |
| *Gallus gallus* | class I | 98.87 % | 98.37 % | 98.89 % | 99.26 % |
| | β2m | 100.00 % | 100.00 % | 100.00 % | 100.00 % |
| *Macaca mulatta* | class I | 99.21 % | 99.04 % | 99.21 % | 100.00 % |
| | β2m | 100.00 % | 100.00 % | 100.00 % | 100.00 % |
| *Rattus norvegicus* | class I | 100.00 % | 100.00 % | 100.00 % | 100.00 % |
| | β2m | 100.00 % | 100.00 % | 100.00 % | 100.00 % |
| *Sus scrofa* | class I | 98.54 % | 98.54 % | 98.54 % | 98.54 % |
| | β2m | 100.00 % | 100.00 % | 100.00 % | 100.00 % |
| **total** | | 99.16 % | 98.92 % | **99.16 %** | **99.39 %** |

Grass carp is the hardest case (one crystal structure; most errors in loops and at the end
of α3). In several of its disputed positions the carp structure and AlphaFold models of
related fish disagree with each other, so part of the remaining difference is uncertainty in
the reference itself.

## MHC: lineages left out

For each lineage, all of its sequences and structures were removed from training, and
AlphaFold DB models of that lineage that pass quality filters (domain pTM ≥ 0.75, mean
pLDDT ≥ 80; residues with pLDDT ≥ 70) were numbered. Reference: the IMGT position each
residue occupies structurally. This reference is derived from structure in a similar way to
the broad set's training data, which favours the broad set somewhat.

| lineage left out | chain | residues | mammal | broad | combined |
|---|---|---|---|---|---|
| ray-finned fish | class I | 155,256 | 91.7 % | 92.9 % | 92.4 % |
| | class II α | 18,631 | 97.2 % | 97.1 % | 97.2 % |
| | class II β | 43,383 | 91.8 % | 93.9 % | 93.6 % |
| | β2m | 4,001 | 98.9 % | 99.0 % | 98.9 % |
| birds | class I | 13,821 | 98.0 % | 98.5 % | 98.0 % |
| | class II α | 2,763 | 95.6 % | 94.2 % | 94.4 % |
| | class II β | 18,635 | 99.1 % | 98.0 % | 98.8 % |
| | β2m | 770 | 99.9 % | 99.7 % | 99.9 % |
| squamates and tuatara | class I | 19,148 | 97.1 % | 97.2 % | 97.1 % |
| | class II α | 1,727 | 96.1 % | 95.9 % | 95.9 % |
| | class II β | 5,936 | 98.7 % | 98.8 % | 98.6 % |
| | β2m | 291 | 99.7 % | 99.7 % | 99.7 % |
| amphibians | class I | 17,586 | 95.1 % | 95.6 % | 95.2 % |
| | class II α | 2,214 | 93.8 % | 94.4 % | 94.3 % |
| | class II β | 4,893 | 98.5 % | 98.3 % | 98.3 % |
| | β2m | 95 | 98.9 % | 98.9 % | 98.9 % |
| mammals | class I | 48,407 | 95.9 % | 95.3 % | 95.8 % |
| | class II α | 4,462 | 97.7 % | 96.7 % | 97.5 % |
| | class II β | 14,046 | 94.2 % | 94.7 % | 94.9 % |
| | β2m | 196 | 98.5 % | 99.5 % | 98.5 % |

Leaving out all mammals also removes human and mouse, the best-studied species, so that row
is a hard stress test rather than a realistic use case.

## TCR constant domains

`imm_tcrnum -const` vs IMGT/3Dstructure-DB numbering of the constant domain, 143 distinct TCR
chains (structures were not used for training; the models come from IMGT/DomainDisplay
sequences). IMGT/3Dstructure-DB numbers the TRAC domain 2001–2128, imm-tools 1001–1128; the
comparison is per IMGT position.

| domain | chains | agreement |
|---|---|---|
| TRAC | 76 | 98.8 % |
| TRBC | 67 | 98.6 % |

The remaining differences are mostly the FG loop, where IMGT/3Dstructure-DB files place
insertions in different ways.

## TCR V domains

V-domain numbering is ANARCI's (IMGT scheme); see Dunbar and Deane, *Bioinformatics* 32,
298–300 (2016).

## Caveats

* Sample sizes differ a lot (e.g. amphibian β2m: 95 residues); one wrongly placed domain can
  move a small row by several percent.
* Errors concentrate in loops and at domain ends, where IMGT itself is least consistent.
* Non-classical MHC (MR1, CD1, HFE, MIC, …), single-chain TCRs and peptides tethered to the
  MHC chain were not specifically tested.
