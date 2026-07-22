#!/bin/bash

export PREFIX="$1"
export BAM="$2"
export VCF="$3"
export REF="$4"
export CONS="$5"

export THRESH_COV=$(bedtools genomecov -bga -ibam "$BAM" | awk '{total+=$4; n++} END {print total*0.05/n}')
#export THRESH_COV=0
bedtools genomecov -bga -ibam "$BAM" | awk '$4 <= var' var="$THRESH_COV" > "$PREFIX.lowcov.bed"
bcftools view "$VCF" -Oz -o "$PREFIX.vcf.gz"
bcftools index "$PREFIX.vcf.gz"
bcftools consensus -f "$REF" -m "$PREFIX.lowcov.bed" "$PREFIX.vcf.gz" > "$CONS"
#bcftools consensus -f "$REF" "$PREFIX.vcf.gz" > "$CONS"