source(snakemake@params$evo_filter)
source(snakemake@params$eval_fun)
source(snakemake@params$get_data)
sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)
#hex_tax <- read.delim(snakemake@input$hex_tax, header=TRUE, sep="\t")
taxonomy <- read.delim(snakemake@input$taxonomy, header=TRUE, sep="\t", row.names=1)
counts <- read.delim(snakemake@input$counts, header=TRUE, sep="\t", row.names=1)
evodistlist <- read.delim(snakemake@input$evodistlist, header=FALSE, sep="\t")
params <- read.delim(snakemake@input$params, header=TRUE, sep="\t")
finbol_tax <- read.delim(snakemake@input$finbol_tax)
se_family <- read.delim(snakemake@input$se_family)
cluster_tax <- get_se_cluster_taxonomy()
i <- snakemake@wildcards$evo_run
ord <- snakemake@wildcards$order

x <- evo_filter(counts,
                evodistlist,
                require_overlap = snakemake@params$require_overlap,
                dist_type       = snakemake@params$dist_type,
                dist_threshold  = snakemake@params$dist_threshold
                )

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
write.table(res, snakemake@output$res_file, sep="\t", row.names=FALSE)
sink()