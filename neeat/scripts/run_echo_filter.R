source(snakemake@params$echo_filter)
source(snakemake@params$eval_fun)
source(snakemake@params$get_data)
sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)
#hex_tax <- read.delim(snakemake@input$hex_tax, header=TRUE, sep="\t")
taxonomy <- read.delim(snakemake@input$taxonomy, header=TRUE, sep="\t", row.names=1)
taxonomy$ASV <- rownames(taxonomy)
counts <- read.delim(snakemake@input$counts, header=TRUE, sep="\t", row.names=1)
matchlist <- read.delim(snakemake@input$matchlist, header=TRUE, sep="\t")
params <- read.delim(snakemake@input$params, header=TRUE, sep="\t")
finbol_tax <- read.delim(snakemake@input$finbol_tax)
se_family <- read.delim(snakemake@input$se_family)
cluster_tax <- get_se_cluster_taxonomy(snakemake@params$cluster_tax)
i <- snakemake@wildcards$i
ord <- snakemake@wildcards$order

x <- echo_filter(counts, 
                 matchlist, 
                 min_match=snakemake@params$min_match,
                 min_overlap=snakemake@params$min_overlap,
                 read_ratio_type=snakemake@params$read_ratio_type,
                 max_read_ratio=snakemake@params$max_read_ratio,
                 require_corr=snakemake@params$require_corr)

outfile <- file(snakemake@output$discarded)
#cat(x$discarded_clusters, sep='\n', file=outfile)
writeLines(x$discarded_clusters, outfile)

y <- eval_res(x$discarded_clusters, taxonomy, finbol_tax, cluster_tax, se_family, debug=FALSE)
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