# SPDX-License-Identifier: Apache-2.0
"""Closest common germline allele of a TCR chain, from pluggable population-frequency sources.

A *source* is one table in models/tcr_alleles/<name>.tsv (built by a script in dev/alleles/).
Every source has the same format, so adding one needs no code here:

    # name: vdjonline                      metadata lines '# key: value'
    # title: VDJ.online gene library       (name, title, species, populations, citation,
    # species: human                        url, licence, built); species: as ANARCI names it
    # populations: all = 134 donors, ...
    gene  allele  segment  population  freq  n  imgt

    gene        IMGT gene name without '/' (TRAV14DV4)
    allele      allele name (TRAV14DV4*02; novel alleles as named by the source)
    segment     V or J
    population  e.g. all, EUR, AFR (a source can give several)
    freq        allele frequency within the gene in that population (0-1)
    n           haplotypes (or individuals) the frequency is based on
    imgt        germline-encoded protein in IMGT numbering: 128 columns, column i = IMGT
                position i, '-' = no residue (V: 1-104, J: 118-128)

The closest common allele: among the alleles of the chain's V (J) gene with freq >= min_freq
in the chosen population, the one with the highest identity over the germline-encoded IMGT
positions of the chain; ties go to the more frequent allele. The CDR3 junction is not
compared (it is trimmed and edited during recombination).
"""
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

DIR = Path(__file__).parent / "models" / "tcr_alleles"
DEFAULT = "vdjonline"                 # source used when none is chosen
DEFAULT_POPULATION = "all"
COMMON_MIN = 0.05                     # "common": frequency >= 5 %
POSITIONS = {"V": range(1, 105), "J": range(118, 129)}
COLUMNS = ("gene", "allele", "segment", "population", "freq", "n", "imgt")


def gene_key(name):
    """'TRAV14/DV4*02' -> 'TRAV14DV4' (gene, without '/' and allele)."""
    return name.split("*")[0].replace("/", "")


@dataclass(frozen=True)
class Allele:
    name: str
    segment: str
    freq: float
    n: int
    imgt: str


@dataclass
class Source:
    name: str
    meta: dict
    alleles: dict = field(default_factory=lambda: defaultdict(list))   # (pop, gene) -> [Allele]

    @property
    def species(self):
        return self.meta.get("species", "human")

    @property
    def populations(self):
        return sorted({p for p, _ in self.alleles})

    def common(self, gene, population, min_freq=COMMON_MIN):
        return [a for a in self.alleles.get((population, gene_key(gene)), ())
                if a.freq >= min_freq]


def available():
    """name -> path of every installed source."""
    return {p.stem: p for p in sorted(DIR.glob("*.tsv"))}


_LOADED = {}


def load(name=DEFAULT):
    """The source `name` (cached). Raises ValueError for an unknown source."""
    if name not in _LOADED:
        path = available().get(name)
        if path is None:
            raise ValueError(f"unknown allele-frequency source '{name}'; "
                             f"available: {', '.join(available()) or 'none'}")
        _LOADED[name] = read(path)
    return _LOADED[name]


def read(path):
    """Parse a source table (format in the module docstring)."""
    meta, src, header = {}, None, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("#"):
            key, _, value = line[1:].partition(":")
            if value:
                meta[key.strip()] = value.strip()
            continue
        if header is None:
            header = line.split("\t")
            if tuple(header) != COLUMNS:
                raise ValueError(f"{path}: columns must be {' '.join(COLUMNS)}")
            src = Source(meta.get("name", Path(path).stem), meta)
            continue
        fields = line.split("\t")
        if len(fields) != len(COLUMNS) or len(fields[-1]) != 128:
            raise ValueError(f"{path}: malformed row: {line[:60]}")
        r = dict(zip(header, fields))
        src.alleles[(r["population"], gene_key(r["gene"]))].append(
            Allele(r["allele"], r["segment"], float(r["freq"]), int(r["n"]), r["imgt"]))
    if src is None:
        raise ValueError(f"{path}: no header line")
    return src


def parse_choice(text):
    """'pmtr' -> ('pmtr', 'all'); 'pmtr/EUR' -> ('pmtr', 'EUR'); '' -> default."""
    name, _, pop = (text or DEFAULT).partition("/")
    return name or DEFAULT, pop or DEFAULT_POPULATION


def check(choice):
    """(source, population) of 'source[/population]'; ValueError if either is unknown."""
    name, population = parse_choice(choice)
    src = load(name)
    if population not in src.populations:
        raise ValueError(f"source '{name}' has no population '{population}'; "
                         f"available: {', '.join(src.populations)}")
    return src, population


def closest(seq, labels, candidates, positions):
    """(Allele, identity %, differences) of the closest candidate, or None.

    seq / labels: the chain and its IMGT numbers ((number, insertion code) or None per residue);
    only `positions` without insertion code are compared (at least 5 shared).
    differences: e.g. ['S57N'] = germline Ser at IMGT 57, Asn in the chain."""
    have = {lab[0]: seq[i] for i, lab in enumerate(labels)
            if lab and not lab[1].strip() and lab[0] in positions}
    best = None
    for a in candidates:
        pairs = [(n, a.imgt[n - 1], have[n]) for n in positions
                 if n in have and a.imgt[n - 1] != "-"]
        if len(pairs) < 5:
            continue
        diff = [f"{g}{n}{x}" for n, g, x in pairs if g != x]
        ident = 100 * (1 - len(diff) / len(pairs))
        if best is None or (ident, a.freq) > (best[1], best[0].freq):
            best = (a, ident, diff)
    return best


EMPTY = dict(closest_common_allele="", common_identity="", common_differences="",
             common_source="")


def annotate(seq, labels, details, choice="", min_freq=COMMON_MIN):
    """The common_* columns of imm_annotate for one TCR chain.

    labels: IMGT numbers of `seq` (from ANARCI); details: ANARCI domain details (v_gene,
    j_gene and their v_species, j_species, else species); choice: 'source' or 'source/population' (default: DEFAULT, 'all')."""
    src, population = check(choice)
    found = []
    for seg in ("V", "J"):
        gene = details.get(f"{seg.lower()}_gene")
        species = details.get(f"{seg.lower()}_species") or details.get("species")
        hit = gene and species == src.species and closest(
            seq, labels, src.common(gene, population, min_freq), POSITIONS[seg])
        if hit:
            found.append(hit)
    if not found:
        return dict(EMPTY)
    return dict(closest_common_allele=" / ".join(a.name for a, _, _ in found),
                common_identity=" / ".join(f"{i:.1f}" for _, i, _ in found),
                common_differences=",".join(d for _, _, ds in found for d in ds) or "none",
                common_source=src.name if population == DEFAULT_POPULATION
                else f"{src.name}/{population}")
