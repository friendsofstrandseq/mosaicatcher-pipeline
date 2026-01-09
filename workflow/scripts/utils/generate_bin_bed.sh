#!/bin/bash
# Generate bin BED file for a reference genome
# Usage: generate_bin_bed.sh <reference_fasta> <output_bed> <window_size>

set -e

REFERENCE_FASTA=$1
OUTPUT_BED=$2
WINDOW_SIZE=${3:-200000}

if [ -z "$REFERENCE_FASTA" ] || [ -z "$OUTPUT_BED" ]; then
    echo "Usage: $0 <reference_fasta> <output_bed> [window_size]"
    echo "Example: $0 canFam3.fa canFam3.bin_200kb_all.bed 200000"
    exit 1
fi

echo "Generating bin BED file for $REFERENCE_FASTA"
echo "Window size: $WINDOW_SIZE bp"
echo "Output: $OUTPUT_BED"

# Index reference if not already indexed
if [ ! -f "${REFERENCE_FASTA}.fai" ]; then
    echo "Indexing reference FASTA..."
    samtools faidx "$REFERENCE_FASTA"
fi

# Extract chromosome sizes
CHROM_SIZES=$(mktemp)
cut -f 1,2 "${REFERENCE_FASTA}.fai" > "$CHROM_SIZES"

# Generate windows
echo "Generating ${WINDOW_SIZE}bp windows..."
bedtools makewindows -g "$CHROM_SIZES" -w "$WINDOW_SIZE" | \
    awk -v OFS="\t" '{print $1, $2, $3, "bin_'$WINDOW_SIZE'kb_"NR}' > "$OUTPUT_BED"

# Cleanup
rm "$CHROM_SIZES"

echo "Generated $(wc -l < "$OUTPUT_BED") bins"
echo "Done: $OUTPUT_BED"
