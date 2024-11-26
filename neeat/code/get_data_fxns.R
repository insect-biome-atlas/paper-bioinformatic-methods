# Functions for getting original data. Modify the paths to fit your system.
get_se_meta <- function(f) { read.delim(f) }

get_se_cluster_taxonomy <- function(f) { read.delim(f) }
get_se_cluster_rep_taxonomy <- function(f) { T<- get_se_cluster_taxonomy(f); T[T$representative==1,] }

get_se_cluster_counts_df <- function(f) { read.delim(f) }
get_se_cluster_counts_dt <- function(f) { require(data.table); fread(f) }

get_se_cluster_taxonomy_samples_hexapoda <- function(meta_f, counts_f, taxonomy_f) {

    meta <- get_se_meta(meta_f)

    samples <- meta$sampleID_NGI[meta$lab_sample_type=="sample" & meta$sequencing_successful==TRUE]

    counts <- get_se_cluster_counts_dt(counts_f)

    index <- match(samples,colnames(counts))
    index <- c(1, index[!is.na(index)])
    counts <- counts[,..index]
    include_rows <- rowSums(counts[,2:ncol(counts)])>0
    counts <- counts[include_rows,]
    include_cols <- c(TRUE,as.logical(colSums(counts[,2:ncol(counts)])>0))
    counts <- counts[,..include_cols]

    taxonomy <- get_se_cluster_taxonomy(taxonomy_f)

    hex_classes <- c("Insecta","Diplura","Protura","Collembola")
    taxonomy <- taxonomy[taxonomy$Class %in% hex_classes,]

    taxonomy[taxonomy$cluster %in% counts$cluster,]
}

