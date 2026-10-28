#!/bin/bash
set -euo pipefail

TMP_DIR="$1"
SEQ_PREFIX="$2"
BAM="$3"
VCF="$4"
REF="$5"
PNG="$6"

TMP_SEQDIR="$TMP_DIR/$(dirname "$SEQ_PREFIX")"
TMP_SEQNAME="$(basename "$SEQ_PREFIX")"
TMP_SEQNAME="${TMP_SEQNAME// /_}"
TMP_REFNAME="$(basename "$REF")"
TMP_REFNAME="${TMP_REFNAME// /_}"

mkdir -p "$TMP_SEQDIR"

safe_cp() {
    src="$1"
    dst="$2"

    if [ "$(realpath "$src")" != "$(realpath -m "$dst")" ]; then
        cp "$src" "$dst"
    fi
}

safe_cp "$BAM" "$TMP_SEQDIR/$TMP_SEQNAME.bam"
safe_cp "$BAM.bai" "$TMP_SEQDIR/$TMP_SEQNAME.bam.bai"
safe_cp "$VCF" "$TMP_SEQDIR/$TMP_SEQNAME.vcf"
safe_cp "$REF" "$TMP_SEQDIR/$TMP_REFNAME"
safe_cp "$REF.fai" "$TMP_SEQDIR/$TMP_REFNAME.fai"

cat <<EOF_BATCH > "${PNG}.txt"
genome $TMP_SEQDIR/$TMP_REFNAME
load $TMP_SEQDIR/$TMP_SEQNAME.bam
load $TMP_SEQDIR/$TMP_SEQNAME.vcf

colorBy READ_STRAND
collapse
squish
setTrackHeight 25
maxPanelHeight 850
snapshotDirectory $TMP_SEQDIR
snapshot $TMP_SEQNAME.igv.png
exit
EOF_BATCH

cat "${PNG}.txt"

xvfb-run \
  --auto-servernum \
  --server-args='-screen 1, 10000x1600x24' \
  java -Duser.home="${TMP_DIR}" \
      -showversion \
      --module-path="/IGV_2.16.2/lib" \
      -Xmx8g @"/IGV_2.16.2/igv.args" \
      --module=org.igv/org.broad.igv.ui.Main "$@" \
      --batch "${PNG}.txt"

safe_cp "$TMP_SEQDIR/$TMP_SEQNAME.igv.png" "$TMP_DIR/$SEQ_PREFIX.igv.png"
safe_cp "$TMP_SEQDIR/$TMP_SEQNAME.igv.png" "$PNG"
