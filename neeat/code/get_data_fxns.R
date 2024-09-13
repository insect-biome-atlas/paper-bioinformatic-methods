# Functions for getting original data. Modify the paths to fit your system.
get_se_meta <- function() { read.delim("~/dev/figshare-repos/iba/raw_data/CO1_sequencing_metadata_SE.tsv", sep="\t") }
get_se_counts <- function() { read.delim("~/dev/figshare-repos/iba/processed_data/SE.v2/cluster_counts.tsv", sep="\t", header=TRUE) }
get_se_cluster_taxonomy <- function() { read.delim("~/dev/figshare-repos/iba/processed_data/SE.v2/cluster_taxonomy.tsv") }


