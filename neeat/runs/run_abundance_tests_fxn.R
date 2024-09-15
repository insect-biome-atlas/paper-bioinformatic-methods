source("../code/abundance_filter.R")
source("../code/evaluate_fxns.R")
source("../code/get_data_fxns.R")

run_abundance_tests <- function(params) {

    data_path <- "../data/"

    hex_tax <- get_se_cluster_taxonomy_samples_hexapoda()
    orders <- unique(hex_tax$Order[!grepl("_X",hex_tax$Order) & !grepl("unclassified",hex_tax$Order)])

    finbol_tax <- read.delim("../evaluation_data/finbol_taxonomy.tsv")
    se_family <- read.delim("../evaluation_data/se_fauna_family.tsv")

    res_file <- paste0("../results/abundance_res_",params$run[1],".tsv")
    res <- data.frame()

    for (ord in orders) {

        cat("Processing",ord,"\n")

        ord_tax <- read.delim(paste0(data_path,ord,"_taxonomy.tsv"), sep="\t", header=TRUE)
        cluster_tax <- hex_tax[hex_tax$Order==ord,]

        # Cycle over runs for this order
        for (i in 1:nrow(params)) {
 
#            print(params[i,])

            if (params$count_type[i] == "raw") {
                ord_counts <- read.delim(paste0(data_path,ord,"_counts.tsv"), sep="\t", header=TRUE, row.names=1)
            } else if (params$count_type[i] == "tot_proportional") {
                ord_counts <- read.delim(paste0(data_path,ord,"_tot_prop_counts.tsv"), sep="\t", header=TRUE, row.names=1)
            } else if (params$count_type[i] == "sample_proportional") {
                ord_counts <- read.delim(paste0(data_path,ord,"_sample_prop_counts.tsv"), sep="\t", header=TRUE, row.names=1)
            } else {
                ord_counts <- read.delim(paste0(data_path,ord,"_cal_counts.tsv"), sep="\t", header=TRUE, row.names=1)
            }

            x <- abundance_filter(ord_counts,
                                  cutoff      = params$cutoff[i],
                                  cutoff_type = params$cutoff_type[i]
                                 )

#            print(str(x))

            y <- eval_res(x$discarded_clusters, ord_tax, finbol_tax, cluster_tax, se_family, debug=FALSE)

#            print(str(y))

            res <- rbind(res,
                         c(taxon=ord,
                           as.list(params[i,]),
                           as.list(y)
                          )
                        )

#            print(str(res))
        }

        # Write results up to this order just in case
        write.table(res,res_file, sep="\t", row.names=FALSE)
    }

    # Compute summary result for all hexapod orders
    res <- add_hexapoda_res(res,params)

    # Write final version of result table
    write.table(res,res_file, sep="\t", row.names=FALSE)
}

