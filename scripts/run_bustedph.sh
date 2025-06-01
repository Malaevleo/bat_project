#!/bin/bash

PAL2NAL_OUTPUT_DIR="batecmreg_sig_pal2nal" # Directory with PAL2NAL alignments of genes undergoing positive selection
PRUNED_TREE_DIR="batecmreg_trees" # Directory with labeled trees for the gene alignments
ABSREL_OUTPUT_DIR="batecmreg_absrel" # Output directory
HYPHY_CMD="hyphy hyphy/tests/hbltests/libv3/BUSTED-PH.bf"

mkdir -p "$ABSREL_OUTPUT_DIR"

for CODON_FILE in "$PAL2NAL_OUTPUT_DIR"/*.fasta; do
    BASENAME=$(basename "$CODON_FILE")
    TREE_FILE="$PRUNED_TREE_DIR/${BASENAME}_pruned.nwk"
    OUTPUT_FILE="$ABSREL_OUTPUT_DIR/${BASENAME}_bustedph.json"
    LOG_FILE="$ABSREL_OUTPUT_DIR/${BASENAME}_bustedph.log"

    if [ -f "$TREE_FILE" ]; then
        # Extract species tags from the alignment file
        ALIGNMENT_SPECIES=$(grep ">" "$CODON_FILE" | cut -c2-3 | sort | uniq)

        # Extract species from the tree
        TREE_SPECIES=$(cat "$TREE_FILE" | grep -oP '[^,():]+' | cut -c1-2 | sort | uniq)

        # Find differences
        DIFF=$(comm -23 <(echo "$ALIGNMENT_SPECIES") <(echo "$TREE_SPECIES"))

        if [ -n "$DIFF" ]; then
            echo "Species mismatch for $CODON_FILE. Missing in tree: $DIFF. Skipping."
            continue
        fi

        echo "Running BUSTED-PH for $CODON_FILE with tree $TREE_FILE..."
        $HYPHY_CMD --alignment "$CODON_FILE" --tree "$TREE_FILE" --output "$OUTPUT_FILE" > "$LOG_FILE" 2>&1

        if [ $? -eq 0 ]; then
            echo "BUSTED-PH analysis completed for $CODON_FILE. Output saved to $OUTPUT_FILE and log saved to $LOG_FILE"
        else
            echo "Error running BUSTED-PH for $CODON_FILE. Check log file: $LOG_FILE"
        fi
    else
        echo "Pruned tree for $CODON_FILE not found. Skipping..."
    fi
done

echo "All BUSTED-PH analyses completed. Results are in $ABSREL_OUTPUT_DIR."
