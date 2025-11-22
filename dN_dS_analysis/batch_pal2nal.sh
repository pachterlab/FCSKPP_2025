#!/bin/bash
#Author : Madhumitha Krishnaswamy

#add in the required paths here

PAL2NAL="" #path to pal2nla.pl 

#note - ensure that cds sequence ids match the alignment sequence ids exactly - otherwise will not work, 
#and ensure no premature stop codons or ambiguous sites in cds file

CDS_DIR="" #directory with cds sequence files
PROT_ALIGN_DIR="" #directory with protein alignment files
OUTPUT_DIR="" #output directory to store the codon alignment in paml format

mkdir -p "$OUTPUT_DIR"
for prot_file in "$PROT_ALIGN_DIR"/*.fasta ; do

    # Skip if no files match the pattern
    [ -e "$prot_file" ] || continue

    filename=$(basename "$prot_file")
    base="${filename%.*}"   
    cds_file="$CDS_DIR/$filename"

    if [[ ! -f "$cds_file" ]]; then
        echo "No CDS file found for: $filename  — skipping"
        continue
    fi

    out_file="$OUTPUT_DIR/${base}.pal2nal"

    echo "Running PAL2NAL on: $filename"

    perl "$PAL2NAL" \
        "$prot_file" \
        "$cds_file" \
        -output paml \
        > "$out_file"

done
