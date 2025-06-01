#!/bin/bash

# Directories
NUC_DIR="batstring_nucs"              # Directory containing nucleotide files
OUTPUT_DIR="batstring_masce"  # Directory to store MACSE output
MACSE_JAR_PATH="macse_v2.07.jar"  # Path to macse.jar (update this)
# Check if MACSE JAR file exists
if [ ! -f "$MACSE_JAR_PATH" ]; then
    echo "Error: MACSE JAR file not found at $MACSE_JAR_PATH"
    exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all .fasta files in the nuc directory
for FILE in "$NUC_DIR"/*.fasta; do
    # Check if the file exists to prevent errors if no files are present
    if [ -f "$FILE" ]; then
        BASENAME=$(basename "$FILE" .fasta)
        OUTPUT_NT="$OUTPUT_DIR/${BASENAME}_aligned_NT.fasta"
        OUTPUT_AA="$OUTPUT_DIR/${BASENAME}_aligned_AA.fasta"

        # Run MACSE with explicit options matching the GUI behavior
        echo "Processing $FILE..."
        java -jar "$MACSE_JAR_PATH" \
            -prog alignSequences \
            -out_AA "$OUTPUT_AA" \
            -out_NT "$OUTPUT_NT" \
            -seq "$FILE"

        # Check if the files were created
        if [ -f "$OUTPUT_NT" ] && [ -f "$OUTPUT_AA" ]; then
            echo "Alignment completed for $FILE."
        else
            echo "Error: No output generated for $FILE. Check MACSE execution."
        fi
    fi
done

echo "All alignments are complete. Check the $OUTPUT_DIR directory."
