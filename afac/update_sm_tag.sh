#!/bin/bash

# Exit immediately if a command exits with a non-zero status,
# treat unset variables as an error, and ensure pipelines fail on the first error.
set -euo pipefail

module load SAMtools/1.14-GCC-11.2.0

SAMPLE_FOLDER=$1

# Define directories
OLD_BAM_DIR="$SAMPLE_FOLDER/old_bam"
NEW_BAM_DIR="$SAMPLE_FOLDER/bam"

# Define the new SM tag value by extracting the sample name from the folder path
NEW_SM_TAG=$(basename "$SAMPLE_FOLDER")

# Create the output directory if it doesn't exist
mkdir -p "$NEW_BAM_DIR"

# Iterate over each BAM file in the old_bam directory
for BAM_FILE in "$OLD_BAM_DIR"/*.bam; do
    echo "Processing $BAM_FILE..."
    # Check if any BAM files exist
    if [[ ! -e "$BAM_FILE" ]]; then
        echo "No BAM files found in $OLD_BAM_DIR. Skipping..."
        break
    fi

    # Extract the base filename (e.g., NA21101_HGSVCpool1NEWiTRU1A68.sort.mdup.bam)
    BASENAME=$(basename "$BAM_FILE")
    echo "Base name: $BASENAME"

    # Capture previous SM tag value
    OLD_SM_TAG=$(samtools view -H "$BAM_FILE" | grep "^@RG" | grep -o "SM:[^[:space:]]*" | cut -d':' -f2)
    echo "Old SM tag: $OLD_SM_TAG"
    echo "New SM tag: $NEW_SM_TAG"

    # Define the output BAM file path
    OUTPUT_BAM="$NEW_BAM_DIR/$BASENAME"
    echo "Output BAM: $OUTPUT_BAM"

    # Create a temporary file for the modified header
    TEMP_HEADER=$(mktemp)

    # Step 1: Extract the header and modify the SM tag
    samtools view -H "$BAM_FILE" |
        sed "s/SM:$OLD_SM_TAG/SM:$NEW_SM_TAG/" >"$TEMP_HEADER"

    # Step 2: Reheader the BAM file with the modified header
    samtools reheader "$TEMP_HEADER" "$BAM_FILE" >"$OUTPUT_BAM"

    # Optional Step: Sort the new BAM file to ensure it is in coordinate order
    # This step can be skipped if the original BAM is already sorted and reheader preserves the order
    samtools sort -o "${OUTPUT_BAM%.bam}.sorted.bam" "$OUTPUT_BAM"
    mv "${OUTPUT_BAM%.bam}.sorted.bam" "$OUTPUT_BAM"

    # Step 3: Index the new BAM file
    samtools index "$OUTPUT_BAM"

    # Clean up the temporary header file
    rm "$TEMP_HEADER"

    echo "Created $OUTPUT_BAM with updated SM tag and indexed it."
done

echo "All BAM files have been processed successfully."
