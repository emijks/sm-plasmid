#!/bin/bash

export TMP_DIR="$1"
export SEQ_PREFIX="$2"
export BAM="$3"
export VCF="$4"
export REF="$5"
export PNG="$6"


export TMP_SEQDIR="$TMP_DIR/$(dirname "$SEQ_PREFIX")"
export TMP_SEQNAME="$(basename "$SEQ_PREFIX")"
export TMP_SEQNAME=${TMP_SEQNAME// /_}
export TMP_REFNAME="$(basename "$REF")"
export TMP_REFNAME=${TMP_REFNAME// /_}


mkdir -p "$TMP_SEQDIR"

cp "$BAM" "$TMP_SEQDIR/$TMP_SEQNAME.bam"
cp "$BAM.bai" "$TMP_SEQDIR/$TMP_SEQNAME.bam.bai"
cp "$VCF" "$TMP_SEQDIR/$TMP_SEQNAME.vcf"
cp "$REF" "$TMP_SEQDIR/$TMP_REFNAME"
cp "$REF.fai" "$TMP_SEQDIR/$TMP_REFNAME.fai"


cat <<EOF > "${PNG}.txt"
genome $TMP_SEQDIR/$TMP_REFNAME
load $TMP_SEQDIR/$TMP_SEQNAME.bam
load $TMP_SEQDIR/$TMP_SEQNAME.vcf
colorBy READ_STRAND
snapshotDirectory $TMP_SEQDIR
maxPanelHeight 3500
snapshot $TMP_SEQNAME.igv.png
exit
EOF

cat "${PNG}.txt"

xvfb-run \
  --auto-servernum \
  --server-args='-screen 1, 7680x4320x24' \
  java -Duser.home="${TMP_DIR}" \
      -showversion \
      --module-path="/IGV_2.16.2/lib" \
      -Xmx8g @"/IGV_2.16.2/igv.args" \
      --module=org.igv/org.broad.igv.ui.Main "$@" \
        --batch "${PNG}.txt"

cp $TMP_SEQDIR/$TMP_SEQNAME.igv.png "$TMP_DIR/$SEQ_PREFIX.igv.png"
cp $TMP_SEQDIR/$TMP_SEQNAME.igv.png "$PNG"
