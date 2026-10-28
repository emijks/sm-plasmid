#!/bin/bash
# igv.sh WORKDIR BAM VCF REF PNG
# Все промежуточные файлы (ссылки на входы, batch-файл, настройки IGV) — в WORKDIR;
# при успешном завершении WORKDIR удаляется.
set -euo pipefail

WORKDIR="$1"
BAM="$2"
VCF="$3"
REF="$4"
PNG="$5"

rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"

# IGV batch не любит пробелы в путях: даём входам простые имена через ссылки
ln -s "$(realpath "$BAM")" "$WORKDIR/sample.bam"
ln -s "$(realpath "$BAM.bai")" "$WORKDIR/sample.bam.bai"
ln -s "$(realpath "$VCF")" "$WORKDIR/sample.vcf"
ln -s "$(realpath "$REF")" "$WORKDIR/reference.fa"
ln -s "$(realpath "$REF.fai")" "$WORKDIR/reference.fa.fai"

cat <<EOF_BATCH > "$WORKDIR/igv_batch.txt"
genome $WORKDIR/reference.fa
load $WORKDIR/sample.bam
load $WORKDIR/sample.vcf

colorBy READ_STRAND
collapse
squish
setTrackHeight 25
maxPanelHeight 850
snapshotDirectory $WORKDIR
snapshot igv.png
exit
EOF_BATCH

cat "$WORKDIR/igv_batch.txt"

xvfb-run \
  --auto-servernum \
  --server-args='-screen 1, 10000x1600x24' \
  java -Duser.home="$WORKDIR" \
      -showversion \
      --module-path="/IGV_2.16.2/lib" \
      -Xmx8g @"/IGV_2.16.2/igv.args" \
      --module=org.igv/org.broad.igv.ui.Main \
      --batch "$WORKDIR/igv_batch.txt"

mv "$WORKDIR/igv.png" "$PNG"
rm -rf "$WORKDIR"
