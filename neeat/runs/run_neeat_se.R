source("../code/neeat_filter.R")
library(data.table)

taxonomy_path <- "~/dev/figshare-repos/iba/processed_data/SE.v2/"
data_path <- "../data/"
res_path <- "neeat_run8_res/"

tax <- read.delim(paste0(taxonomy_path,"cluster_taxonomy.tsv"), sep="\t", header=TRUE)
tax <- tax[tax$representative==1,]
counts <- fread("../bigdata/cluster_counts.tsv")
orders <- unique(tax$Order[!grepl("_X",tax$Order) & !grepl("unclassified",tax$Order) & (tax$cluster %in% counts$cluster)])
rm(counts)

for (ord in orders) {

    cat("Processing",ord,"\n")

    ord_counts <- read.delim(paste0(data_path,ord,"_cal_counts.tsv"), sep="\t", header=TRUE, row.names=1)
    ord_matchlist <- read.delim(paste0(data_path,ord,"_evodistlist.tsv"))
    ord_tax <- read.delim(paste0(data_path,ord,"_taxonomy.tsv"), sep="\t", header=TRUE)

    res <- neeat_filter(ord_counts, ord_matchlist, ord_tax, debug=TRUE, dump_prefix=paste0(res_path,"All_orders"))

    if (length(res$discarded_clusters)>0)
        cat(res$discarded_clusters,file=paste0(res_path,"neeat_discarded_clusters_se.tsv"),sep="\n",append=TRUE)
    if (length(res$retained_clusters)>0)
        cat(res$retained_clusters,file=paste0(res_path,"neeat_retained_clusters_se.tsv"),sep="\n",append=TRUE)
    if (ord==orders[1])
        write.table(res$filtered_counts,file=paste0(res_path,"filtered_counts_se.tsv"),sep="\t",row.names=FALSE)
    else if (length(res$retained_clusters)>0)
        write.table(res$filtered_counts,file=paste0(res_path,"filtered_counts_se.tsv"),sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE)
}

