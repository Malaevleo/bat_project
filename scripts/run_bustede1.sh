#!/bin/bash

PAL2NAL_OUTPUT_DIR="batstring_pal2nal" # Directory with PAL2NAL alignments of genes
PRUNED_TREE_DIR="batstring_trees" # Directory with labeled trees for gene alignments
ABSREL_OUTPUT_DIR="batstring_bustede1" # Output directory
HYPHY_CMD="hyphy busted"

mkdir -p "$ABSREL_OUTPUT_DIR"

for CODON_FILE in "$PAL2NAL_OUTPUT_DIR"/*.fasta; do
    BASENAME=$(basename "$CODON_FILE")
    TREE_FILE="$PRUNED_TREE_DIR/${BASENAME}_pruned.nwk"
    OUTPUT_FILE="$ABSREL_OUTPUT_DIR/${BASENAME}_busted.json"
    LOG_FILE="$ABSREL_OUTPUT_DIR/${BASENAME}_busted.log"

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

        echo "Running BUSTED for $CODON_FILE with tree $TREE_FILE..."
        $HYPHY_CMD --alignment "$CODON_FILE" --tree "$TREE_FILE" --branches Foreground --error-sink Yes --output "$OUTPUT_FILE" > "$LOG_FILE" 2>&1

        if [ $? -eq 0 ]; then
            echo "BUSTED analysis completed for $CODON_FILE. Output saved to $OUTPUT_FILE and log saved to $LOG_FILE"
        else
            echo "Error running BUSTED for $CODON_FILE. Check log file: $LOG_FILE"
        fi
    else
        echo "Pruned tree for $CODON_FILE not found. Skipping..."
    fi
done

echo "All BUSTED analyses completed. Results are in $ABSREL_OUTPUT_DIR."
