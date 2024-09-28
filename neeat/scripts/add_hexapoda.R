source(snakemake@params$eval_fun)
source(snakemake@params$get_data)
sink(file = snakemake@log[[1]], append = FALSE, type = c("output", "message"),
     split = FALSE)
hex_tax <- read.delim(snakemake@input$hex_tax, header=TRUE, sep="\t")
params <- read.delim(snakemake@input$params, header=TRUE, sep="\t")
finbol_tax <- read.delim(snakemake@input$finbol_tax)
se_family <- read.delim(snakemake@input$se_family)
cluster_tax <- get_se_cluster_taxonomy()
res <- read.delim(snakemake@input$res, header=TRUE, sep="\t")

res <- add_hexapoda_res(res, params)
write.table(res, snakemake@output$res_file, sep="\t", row.names=FALSE)
sink()