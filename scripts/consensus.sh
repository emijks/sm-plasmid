#!/bin/bash
# Использование: consensus.sh WORKDIR BAM VCF REF CONS
# Все промежуточные файлы хранятся в WORKDIR.
# При успешном завершении WORKDIR удаляется.

set -euo pipefail

if [[ $# -ne 5 ]]; then
    printf 'Использование: %s WORKDIR BAM VCF REF CONS\n' "$0" >&2
    exit 2
fi

WORKDIR="$1"
BAM="$2"
VCF="$3"
REF="$4"
CONS="$5"

rm -rf -- "$WORKDIR"
mkdir -p -- "$WORKDIR"

bedtools genomecov -bga -ibam "$BAM" > "$WORKDIR/cov.bedgraph"

if [[ ! -s "$WORKDIR/cov.bedgraph" ]]; then
    printf 'Ошибка: bedtools genomecov вернул пустой файл для %s\n' "$BAM" >&2
    exit 1
fi

THRESH_COV=$(awk '{ total += $4 } END { print total * 0.05 / NR }' \
    "$WORKDIR/cov.bedgraph")

awk -v threshold="$THRESH_COV" '
    BEGIN { FS = OFS = "\t" }
    $4 <= threshold { print $1, $2, $3 }
' "$WORKDIR/cov.bedgraph" > "$WORKDIR/lowcov_raw.bed"

# Исключаем из маски REF-интервалы вызванных инделов.
# Удалённые основания genomecov может считать непокрытыми.
awk '
    BEGIN { FS = OFS = "\t" }
    !/^#/ && length($4) != length($5) {
        print $1, $2 - 1, $2 - 1 + length($4)
    }
' "$VCF" > "$WORKDIR/indels.bed"

if [[ ! -s "$WORKDIR/lowcov_raw.bed" || ! -s "$WORKDIR/indels.bed" ]]; then
    cp -- "$WORKDIR/lowcov_raw.bed" "$WORKDIR/lowcov.bed"
else
    # Один словарь для обоих BED-файлов.
    LC_ALL=C awk \
        -v low="$WORKDIR/lowcov.ascii.bed" \
        -v indels="$WORKDIR/indels.ascii.bed" '
        BEGIN { FS = OFS = "\t" }
        !($1 in id) {
            id[$1] = "ctg" ++n
            print id[$1], $1
        }
        {
            output = (FILENAME == ARGV[1] ? low : indels)
            print id[$1], $2, $3 > output
        }
    ' "$WORKDIR/lowcov_raw.bed" "$WORKDIR/indels.bed" \
        > "$WORKDIR/contigs.tsv"

    bedtools subtract \
        -a "$WORKDIR/lowcov.ascii.bed" \
        -b "$WORKDIR/indels.ascii.bed" |
        LC_ALL=C awk '
            BEGIN { FS = OFS = "\t" }
            FILENAME == ARGV[1] { name[$1] = $2; next }
            { print name[$1], $2, $3 }
        ' "$WORKDIR/contigs.tsv" - > "$WORKDIR/lowcov.bed"
fi

bcftools view "$VCF" -Oz -o "$WORKDIR/variants.vcf.gz"
bcftools index "$WORKDIR/variants.vcf.gz"

# Маска задана в координатах референса и применяется внутри consensus.
consensus_args=(-f "$REF")
if [[ -s "$WORKDIR/lowcov.bed" ]]; then
    consensus_args+=(-m "$WORKDIR/lowcov.bed")
fi

bcftools consensus "${consensus_args[@]}" "$WORKDIR/variants.vcf.gz" \
    > "$WORKDIR/consensus.fa"

mv -- "$WORKDIR/consensus.fa" "$CONS"
rm -rf -- "$WORKDIR"