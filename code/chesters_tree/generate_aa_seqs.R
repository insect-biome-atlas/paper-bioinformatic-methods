# Generate unaligned aa sequence files

# Set paths
data_path <- "../../data/chesters_tree/"
infile <- "chesters_new_outgroups.fasta"
outfile <- "chesters_new_outgroups_aa.fasta"

# Libraries needed
library(ape)

# Read in fasta sequences
seqs <- read.FASTA(paste0(data_path,infile))

# Note: the trans function does not quite handle ambiguities correctly, generating a few
# warnings in the pal2nal step. This should have negligible effect on the analysis.
aa_seqs <- trans(seqs, code=5)

write.FASTA(aa_seqs, file=paste0(data_path,outfile))

