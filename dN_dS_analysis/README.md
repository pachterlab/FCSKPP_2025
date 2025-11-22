Author: Madhumitha Krishnaswamy

This contains all the relevant code and data required to reproduce supplementary figure S8 and perform dN/dS calculations.

Data:
1. relevant_gene_expression_sorted.csv - list of 167 genes and expression data, sorted on the basis of median expression level
2. all_genes.txt  - list of 167 human gene ids (Ensembl) used for analysis - subset from data file above
3. gene_expression_dnds_159_genes.csv - final data file used for generating figure

Codes:
Run the below codes in order (for all the shell scripts just run bash script.sh), I recommend keeping the file ids the same as the gene ids for easier use downstream, and save all the outputs in different directories 

1. ensembl_sequence_fetch.ipynb - Fetching the sequence files using ensembl rest api for the given genes - using all_genes.txt, and filter on the basis of species, taking only single copy representative. Save the corresponding cds and protein sequences in separate subdirectories, with the same sequence ID header and file name in all.
2. batch_mafft.sh  - run on directory with protein sequences to get the alignment files - can modify alignment parameters depending on number of sequences, recommend checking for gaps as gappy regions can be problematic during codon alignment and dn/ds calculations.
3. batch_pal2nal.sh - performs codon alignment and saves in format compatible for paml analysis - requires alignment and cds as input.
4. batch_codeml.sh - dN/dS calculations, with the control file written in the script - can be modified for different calculations, saves a summary txt.
5. plot_dnds_R.ipynb - uses the final output file to generate boxplot using ggplot2 with significance values using ggpubr.