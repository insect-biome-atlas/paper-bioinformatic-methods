#! /bin/bash
# Construct aligned data matrix in Nexus format for
# Entognatha and outgroup analysis with MrBayes

# Merge taxonomies
R --no-save < merge_taxonomies.R

# Read coding DNA sequences and get aa sequences
# The R script has hardcoded infile path = "../../data/entognatha_outgroup/entognatha_outgroup.fasta"
# and outfile path = "../../data/entognatha_outgroup/entognatha_outgroup_aa.fasta"
# Note some minor errors in the handling of ambiguities ("Ns")
R --no-save < generate_aa_seqs.R

# Align aa sequences using mafft --auto
mafft --auto ../../data/entognatha_outgroup/entognatha_outgroup_aa.fasta > ../../data/entognatha_outgroup/entognatha_outgroup_aa_aligned.fasta

# Convert them back to nucleotide sequences
# pal2nal will correctly report a few errors in the handling of ambiguities. This should have negligible effect on
# the analysis and is ignored here.
pal2nal.pl ../../data/entognatha_outgroup/entognatha_outgroup_aa_aligned.fasta ../../data/entognatha_outgroup/entognatha_outgroup.fasta -codontable 5 -output fasta > ../../data/entognatha_outgroup/entognatha_outgroup_aligned.fasta

# Make the NEXUS constraint files
R --no-save < generate_constraint_files.R

