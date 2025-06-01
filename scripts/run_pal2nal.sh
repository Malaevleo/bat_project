#!/bin/bash

# Directories
INPUT_DIR="batstring_masce"          # Directory containing aligned AA and NT files
OUTPUT_DIR="batstring_pal2nal"       # Directory to store PAL2NAL results
PAL2NAL_PATH="pal2nal.v14/pal2nal.pl" # Path to the PAL2NAL script (update this)

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all amino acid files in the input directory
for AA_FILE in "$INPUT_DIR"/*_AA.fasta; do
    # Derive the corresponding nucleotide file
    BASENAME=$(basename "$AA_FILE" "_aligned_AA.fasta")
    NT_FILE="$INPUT_DIR/${BASENAME}_aligned_NT.fasta"

    # Check if the nucleotide file exists
    if [ -f "$NT_FILE" ]; then
        # Define the output file
        OUTPUT_FILE="$OUTPUT_DIR/${BASENAME}_codon_aligned.fasta"

        # Run PAL2NAL
        echo "Running PAL2NAL for $BASENAME..."
        perl "$PAL2NAL_PATH" "$AA_FILE" "$NT_FILE" -output fasta -nogap > "$OUTPUT_FILE"

        if [ $? -eq 0 ]; then
            echo "Codon alignment completed: $OUTPUT_FILE"
        else
            echo "Error running PAL2NAL for $BASENAME."
        fi
    else
        echo "Nucleotide file not found for $BASENAME. Skipping..."
    fi
done

echo "All PAL2NAL runs completed. Results are in $OUTPUT_DIR."
