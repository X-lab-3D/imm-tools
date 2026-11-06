#!/usr/bin/env bash
# Standardise a folder of TCR-pMHC structures: first complex only, standard chain IDs,
# IMGT numbering everywhere; then one annotation table and a list of anything imm_check
# does not like.
#
#   bash batch.sh structures/ out/
#
# out/<name>_imgt.pdb (or .cif)   standardised structures
# out/annotation.tsv              imm_annotate of all of them (one header line)
# out/problems.tsv                imm_check rows that are not OK (empty = all fine)
# out/failed.txt                  files that could not be processed, with the reason
set -u -o pipefail
in=${1:?usage: batch.sh <input folder> <output folder>}
out=${2:?usage: batch.sh <input folder> <output folder>}
mkdir -p "$out"
: > "$out/failed.txt"

for f in "$in"/*.pdb "$in"/*.cif; do
  [ -e "$f" ] || continue
  name=$(basename "${f%.*}"); ext=${f##*.}
  if ! { imm_complex -1 "$f" | imm_rechain | imm_tcrnum -const | imm_mhcnum | imm_pepnum \
         > "$out/${name}_imgt.$ext"; } 2> "$out/.err"; then
    echo -e "$name\t$(grep -m1 ERROR "$out/.err")" >> "$out/failed.txt"
    rm -f "$out/${name}_imgt.$ext"
  fi
done
rm -f "$out/.err"

first=1
for f in "$out"/*_imgt.*; do
  [ -e "$f" ] || continue
  if [ $first = 1 ]; then imm_annotate "$f"; first=0; else imm_annotate "$f" | tail -n +2; fi
done > "$out/annotation.tsv" 2>/dev/null

for f in "$out"/*_imgt.*; do
  [ -e "$f" ] || continue
  imm_check "$f" 2>/dev/null | awk -F'\t' 'NR > 1 && $5 != "OK"'
done > "$out/problems.tsv"

echo "$(ls "$out"/*_imgt.* 2>/dev/null | wc -l) structures in $out;" \
     "$(wc -l < "$out/problems.tsv") check warnings; $(wc -l < "$out/failed.txt") failed" >&2
