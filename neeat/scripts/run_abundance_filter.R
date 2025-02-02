source(snakemake@params$abundance_filter)
source(snakemake@params$eval_fun)
source(snakemake@params$get_data)
sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)

taxonomy <- read.delim(snakemake@input$taxonomy, header=TRUE, sep="\t", row.names=1)
taxonomy$ASV <- rownames(taxonomy)
counts <- read.delim(snakemake@input$counts, header=TRUE, sep="\t", row.names=1)
params <- read.delim(snakemake@input$params, header=TRUE, sep="\t")
finbol_tax <- read.delim(snakemake@input$finbol_tax)
se_family <- read.delim(snakemake@input$se_family)
cluster_tax <- get_se_cluster_taxonomy(snakemake@params$cluster_tax)
i <- snakemake@wildcards$abundance_run
ord <- snakemake@wildcards$order

x <- abundance_filter(counts,
                      cutoff = snakemake@params$cutoff,
                      cutoff_type = snakemake@params$cutoff_type
                      )

discarded_otus <- x$discarded_clusters
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