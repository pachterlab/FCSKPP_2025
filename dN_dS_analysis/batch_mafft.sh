#!/bin/bash
#Author : Madhumitha Krishnaswamy

#add in the required paths here

INPUT_DIR=""        # directory with species FASTA files
OUTPUT_DIR=""        # directory to write alignments
MAFFT_BIN="/programs/mafft/bin/mafft"   # path for mafft exe

mkdir -p "$OUTPUT_DIR"

for fasta in "$INPUT_DIR"/*.fasta; do
    fname=$(basename "$fasta")
    aligned_file="$OUTPUT_DIR/$fname"
    "$MAFFT_BIN" "$fasta" > "$aligned_file"
done