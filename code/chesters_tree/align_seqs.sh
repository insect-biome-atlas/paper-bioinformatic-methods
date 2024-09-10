#! /bin/bash
# Construct aligned data matrix for phylogenetic annotation
# with the Chesters tree with new outgroups

# Set data path
data_path="../../data/chesters_tree/"

# Read coding DNA sequences and translate to aa sequences
R --no-save < generate_aa_seqs.R

# Align aa sequences using mafft --auto
mafft --auto ${data_path}chesters_new_outgroups_aa.fasta > ${data_path}chesters_new_outgroups_aa_aligned.fasta

# Convert them back to nucleotide sequences
# pal2nal will correctly report a few errors in the handling of ambiguities. This should have negligible effect on
# the analysis and is ignored here.
pal2nal.pl ${data_path}chesters_new_outgroups_aa_aligned.fasta ${data_path}chesters_new_outgroups.fasta -codontable 5 -output fasta > ${data_path}chesters_new_outgroups_aligned.fasta

gzip ${data_path}chesters_new_outgroups_aligned.fasta

