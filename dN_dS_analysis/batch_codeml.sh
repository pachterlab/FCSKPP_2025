#!/bin/bash
#Author : Madhumitha Krishnaswamy

#add in the required paths here

CODON_DIR=""     # folder with pal2nal codon alignments
TREE_FILE=""       # reference tree for all genes
OUTPUT_DIR=""      # folder to save codeml outputs
SUMMARY_FILE="" # summary table
CODEML_EXE=""      # path to codeml executable

mkdir -p "$OUTPUT_DIR"

echo -e "Gene\tOmega\tlnL" > "$SUMMARY_FILE"

for aln_file in "$CODON_DIR"/*.pal2nal; do
    [ -e "$aln_file" ] || continue 

    gene_name=$(basename "$aln_file" .pal2nal)
    out_file="$OUTPUT_DIR/${gene_name}.out"
    ctl_file="$OUTPUT_DIR/${gene_name}.ctl"

    # Codeml control file - modify based on parameters to be used
    cat > "$ctl_file" << EOF
seqfile = $aln_file
treefile = $TREE_FILE
outfile = $out_file
runmode = 0
seqtype = 1
model = 0
NSsites = 0
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
cleandata = 1
clock = 1
method = 1
EOF

    echo "Running codeml for $gene_name..."
    "$CODEML_EXE" "$ctl_file"


    omega=$(grep "omega (dN/dS) =" "$out_file" | awk -F= '{print $2}' | tr -d ' ')
    lnL=$(grep "^lnL(" "$out_file" | awk '{print $5}')

    
    echo -e "${gene_name}\t${omega}\t${lnL}" >> "$SUMMARY_FILE"

done