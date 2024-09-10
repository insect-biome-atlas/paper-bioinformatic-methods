# Generate unaligned aa sequence files

data_path <- "../../data/entognatha_outgroup/"
infile <- "entognatha_outgroup.fasta"
outfile <- "entognatha_outgroup_aa.fasta"

# Libraries needed
library(ape)

# Read in fasta sequences
cat("Reading nuc sequences\n")
seqs <- read.FASTA(paste0(data_path,infile))

# Note: the trans function does not quite handle ambiguities correctly, generating a few
# warnings in the pal2nal step. This should have negligible effect on the analysis.
aa_seqs <- trans(seqs, code=5)

cat("Writing aa sequences\n")
write.FASTA(aa_seqs, file=paste0(data_path,outfile))


