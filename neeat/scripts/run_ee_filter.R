source(snakemake@params$evo_filter)
source(snakemake@params$eval_fun)
source(snakemake@params$get_data)
sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

taxonomy <- read.delim(snakemake@input$taxonomy, header=TRUE, sep="\t", row.names=1)
taxonomy$ASV <- rownames(taxonomy)
counts <- read.delim(snakemake@input$counts, header=TRUE, sep="\t", row.names=1)
evodistlist <- read.delim(snakemake@input$evodistlist, header=TRUE, sep="\t")
params <- read.delim(snakemake@input$params, header=TRUE, sep="\t")
finbol_tax <- read.delim(snakemake@input$finbol_tax)
se_family <- read.delim(snakemake@input$se_family)
cluster_tax <- get_se_cluster_taxonomy(snakemake@params$cluster_tax)
i <- snakemake@wildcards$ee_run
ord <- snakemake@wildcards$order

if (file.size(snakemake@input$discarded) > 0) {
     echo_set <- read.delim(snakemake@input$discarded, header=FALSE)
     echo_set <- echo_set$V1[echo_set$V1 %in% taxonomy$ASV]
     counts_subset <- counts[!(rownames(counts) %in% echo_set),]
     evodistlist_subset <- evodistlist[!(evodistlist$asv1 %in% echo_set) & !(evodistlist$asv2 %in% echo_set),]
} else {
     echo_set <- character(0)
     counts_subset <- counts
     evodistlist_subset <- evodistlist
}


x <- evo_filter(counts_subset,
                evodistlist_subset,
                require_overlap = snakemake@params$require_overlap,
                dist_type = snakemake@params$dist_type,
                dist_threshold = snakemake@params$dist_threshold
                           )
discarded_otus <- c(echo_set, x$discarded_clusters)
outfile <- file(snakemake@output$discarded)
#cat(x$discarded_clusters, sep='\n', file=outfile)
writeLines(discarded_otus, outfile)

y <- eval_res(discarded_otus, taxonomy, finbol_tax, cluster_tax, se_family, debug=FALSE)

res <- as.data.frame(
     c(taxon=ord, 
          as.list(params[i,]),
          as.list(y)
     )
     )
# Do this after all parameter runs have been completed
#res <- add_hexapoda_res(res, params[i,])
write.table(res, snakemake@output$res_file, sep="\t", row.names=FALSE)
sink()