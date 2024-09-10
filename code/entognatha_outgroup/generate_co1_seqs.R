library(ape)
source("../seq_fxns.R")

# Set data path
data_path <- "../../data/entognatha_outgroup/"


# Generate root sequences
# -----------------------

seqs <- read.FASTA(paste0(data_path,"root_mt_genomes.fasta"))
meta <- read.delim(paste0(data_path,"root_taxa.csv"),sep=";")

# Relabel sequences with species names
root_seqs <- rename_seqs(seqs, meta)

# Extract the coding part
root_seqs <- extract_coding(root_seqs, meta)


# Generate extra entognatha sequences
# -----------------------------------

seqs <- read.FASTA(paste0(data_path,"co1_missing_entognatha_seqs.fasta"))
meta <- read.delim(paste0(data_path,"entognatha_taxa.csv"),sep=";")

# Relabel sequences with species names
ento_seqs <- rename_seqs(seqs, meta)

# Extract the coding part
ento_seqs <- extract_coding(ento_seqs, meta)


# Read in Collembola sequences
seqs <- read.FASTA(paste0(data_path,"co1_collembola_seqs.fasta"))
meta <- read.delim(paste0(data_path,"collembola_taxa.csv"),sep=";")

# Sequences should already be named correctly with tip labels

# Extract the coding part
coll_seqs <- extract_coding(seqs, meta)


# Write the resulting sequences
# -----------------------------

out_file <- paste0(data_path,"entognatha_outgroup.fasta")
write.FASTA(coll_seqs, out_file)
write.FASTA(ento_seqs, out_file, append=TRUE)
write.FASTA(root_seqs, out_file, append=TRUE)


