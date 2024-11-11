# Generate dataset partitions
library(data.table)
source("/cfs/klemming/projects/snic/snic2020-16-248/git/paper-bioinformatic-methods/neeat/code/get_data_fxns.R")

output_path <- "/cfs/klemming/projects/snic/snic2020-16-248/git/paper-bioinformatic-methods/neeat/data/"

taxonomy <- get_se_cluster_rep_taxonomy()
counts <- fread("/cfs/klemming/projects/snic/snic2020-16-248/git/paper-bioinformatic-methods/neeat/bigdata/cluster_counts.tsv", header=TRUE, sep="\t")
cal_counts <- fread("/cfs/klemming/projects/snic/snic2020-16-248/git/paper-bioinformatic-methods/neeat/bigdata/calibrated_cluster_counts.tsv", header=TRUE, sep="\t")
tprop_counts <- fread("/cfs/klemming/projects/snic/snic2020-16-248/git/paper-bioinformatic-methods/neeat/bigdata/tot_proportional_cluster_counts.tsv", header=TRUE, sep="\t")
sprop_counts <- fread("/cfs/klemming/projects/snic/snic2020-16-248/git/paper-bioinformatic-methods/neeat/bigdata/sample_proportional_cluster_counts.tsv", header=TRUE, sep="\t")
#matchlist <- read.delim("neeat/bigdata/vsearch_matchlist.tsv",header=FALSE,sep="\t")

orders <- unique(taxonomy$Order[!grepl("_X",taxonomy$Order) & !grepl("unclassified",taxonomy$Order) & taxonomy$cluster %in% counts$cluster])

for (ord in orders) {
    cat("Processing order ", ord, "\n")
    ord_tax <- taxonomy[taxonomy$Order==ord,]
    # Skip orders with only one cluster
    if (nrow(ord_tax) < 2) {
        cat("Skipping order ", ord, " with only one cluster\n")
        next
    }
    ord_counts <- counts[counts$cluster %in% ord_tax$cluster,]
    ord_counts$cluster <- ord_tax$ASV[match(ord_counts$cluster,ord_tax$cluster)]
    file_name <- paste0(output_path,ord,"_counts.tsv")
    write.table(ord_counts, file_name, row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

    # Update ord_tax to only those clusters that have counts (pos controls are in taxonomy, generating extra names)
    # Note that cluster name is now ASV id in counts
    ord_tax <- ord_tax[ord_tax$ASV %in% ord_counts$cluster,]
    file_name <- paste0(output_path,ord,"_taxonomy.tsv")
    write.table(ord_tax, file_name, sep="\t",row.names=FALSE)

    ord_counts <- cal_counts[cal_counts$cluster %in% ord_tax$cluster,]
    ord_counts$cluster <- ord_tax$ASV[match(ord_counts$cluster,ord_tax$cluster)]
    file_name <- paste0(output_path,ord,"_cal_counts.tsv")
    write.table(ord_counts, file_name, row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

    ord_counts <- tprop_counts[tprop_counts$cluster %in% ord_tax$cluster,]
    ord_counts$cluster <- ord_tax$ASV[match(ord_counts$cluster,ord_tax$cluster)]
    file_name <- paste0(output_path,ord,"_tot_prop_counts.tsv")
    write.table(ord_counts, file_name, row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

    ord_counts <- sprop_counts[sprop_counts$cluster %in% ord_tax$cluster,]
    ord_counts$cluster <- ord_tax$ASV[match(ord_counts$cluster,ord_tax$cluster)]
    file_name <- paste0(output_path,ord,"_sample_prop_counts.tsv")
    write.table(ord_counts, file_name, row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

    #ord_matchlist <- matchlist[(matchlist$V1 %in% ord_tax$cluster) & (matchlist$V2 %in% ord_tax$cluster),]
    #file_name <- paste0(output_path,ord,"_evodistlist.tsv")
    #write.table(ord_matchlist, file_name, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
}

