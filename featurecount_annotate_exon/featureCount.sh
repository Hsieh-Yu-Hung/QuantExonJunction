#!/bin/bash

BAM_FILE=$1
GENOME_GTF=$2
OUTPUT_PREFIX=$3

featureCounts \
 -T 8 \
 -a $GENOME_GTF \
 -o ${OUTPUT_PREFIX} \
 -t exon -g transcript_id \
 -M -J -p -C --countReadPairs \
 -B $BAM_FILE

# Convert FeatureCount Result .jcounts file to bed
echo "Convert FeatureCount Result .jcounts file to bed"
awk -v OFS="\t" '$3 == $6 {print $3,$4,$7,$1,$9}' ${OUTPUT_PREFIX}.jcounts > ${OUTPUT_PREFIX}.junction.bed