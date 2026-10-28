#!/bin/bash
# consensus.sh WORKDIR BAM VCF REF CONS
# Все промежуточные файлы — в WORKDIR; при успешном завершении WORKDIR удаляется.

set -euo pipefail

WORKDIR="$1"
BAM="$2"
VCF="$3"
REF="$4"
CONS="$5"

rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"

bedtools genomecov -bga -ibam "$BAM" > "$WORKDIR/cov.bedgraph"
THRESH_COV=$(awk '{total+=$4; n++} END {print total*0.05/n}' "$WORKDIR/cov.bedgraph")
awk '$4 <= var' var="$THRESH_COV" "$WORKDIR/cov.bedgraph" > "$WORKDIR/lowcov_raw.bed"

# Удалённые основания genomecov не считает покрытыми, поэтому интервал
# настоящей делеции всегда попадает в маску, и bcftools consensus тогда
# пропустил бы делецию и поставил N. Исключаем из маски REF-интервалы
# вызванных инделов: делеция подтверждена ридами (фильтры caller).
awk 'BEGIN {OFS="\t"} !/^#/ && length($4) != length($5) {print $1, $2 - 1, $2 - 1 + length($4)}' \
    "$VCF" > "$WORKDIR/indels.bed"
bedtools subtract -a "$WORKDIR/lowcov_raw.bed" -b "$WORKDIR/indels.bed" > "$WORKDIR/lowcov.bed"

bcftools view "$VCF" -Oz -o "$WORKDIR/variants.vcf.gz"
bcftools index "$WORKDIR/variants.vcf.gz"

# Маска задаётся в координатах референса, поэтому применяется внутри
# bcftools consensus (-m), а не после него: после инделов координаты
# консенсуса сдвигаются относительно референса.
bcftools consensus -f "$REF" -m "$WORKDIR/lowcov.bed" "$WORKDIR/variants.vcf.gz" > "$WORKDIR/consensus.fa"

mv "$WORKDIR/consensus.fa" "$CONS"
rm -rf "$WORKDIR"
