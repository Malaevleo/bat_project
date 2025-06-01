#!/bin/bash

SOURCE_DIR="batecmreg_pal2nal"                # Directory with all alignment files
DEST_DIR="batecmreg_sig_pal2nal"               # Destination directory for significant alignments
SIGNIFICANT_GENES_FILE="batecmreg_bustede1_FDR.txt"  # File with significant gene IDs

mkdir -p "$DEST_DIR"

echo "=== Printing raw contents of significant genes file ==="
while IFS= read -r line; do
    echo "[$line]"
done < "$SIGNIFICANT_GENES_FILE"

# Read significant gene IDs, remove carriage returns and extra whitespace
mapfile -t sig_genes < <(sed 's/\r//g; s/^[ \t]*//; s/[ \t]*$//' "$SIGNIFICANT_GENES_FILE")

echo "=== Processed significant genes (after trimming) ==="
for gene in "${sig_genes[@]}"; do
    echo "[$gene]"
done

# Loop through all alignment files
for alignment in "$SOURCE_DIR"/*aligned_codon_alignment.fasta; do
    base=$(basename "$alignment")
    
    # Extract the gene ID (assumes it's everything before the first underscore)
    gene_id=$(echo "$base" | cut -d'_' -f1-2 | xargs)
    echo "----"
    echo "Processing file: [$base]"
    echo "Extracted gene ID: [$gene_id]"
    
    match_found=false
    for sig in "${sig_genes[@]}"; do
        # Remove any potential carriage returns from the significant gene string
        sig_clean=$(echo "$sig" | tr -d '\r' | xargs)
        echo "Comparing with significant gene: [$sig_clean]"
        if [ "$gene_id" == "$sig_clean" ]; then
            match_found=true
            break
        fi
    done
    
    if [ "$match_found" = true ]; then
        echo "Copying $base because gene ID [$gene_id] is significant."
        cp "$alignment" "$DEST_DIR"
    else
        echo "Skipping $base (gene ID [$gene_id] is not significant)."
    fi
done

echo "Finished copying significant gene alignments to $DEST_DIR."
