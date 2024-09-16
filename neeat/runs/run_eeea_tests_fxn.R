source("../code/abundance_filter.R")
source("../code/evaluate_fxns.R")
source("../code/get_data_fxns.R")

run_eeea_tests <- function(params) {

    data_path <- "../data/"
    eval_data_path <- "../evaluation_data/"
    raw_res_path <- "../raw_results/"

    hex_tax <- get_se_cluster_taxonomy_samples_hexapoda()
    orders <- unique(hex_tax$Order[!grepl("_X",hex_tax$Order) & !grepl("unclassified",hex_tax$Order)])

    se_family <- read.delim(paste0(eval_data_path,"se_fauna_family.tsv"))
    finbol_tax <- read.delim(paste0(eval_data_path,"finbol_taxonomy.tsv"))
    cluster_tax <- get_se_cluster_taxonomy()

    res_file <- paste0("../results/eeea_res_",params$run[1],"-",params$run[nrow(params)],".tsv")
    res <- data.frame()

    for (ord in orders) {

        cat("Processing",ord,"\n")

        ord_tax <- read.delim(paste0(data_path,ord,"_taxonomy.tsv"), sep="\t", header=TRUE)

        ord_counts <- read.delim(paste0(data_path,ord,"_cal_counts.tsv"), sep="\t", header=TRUE, row.names=1)

        # Cycle over runs for this order
        for (i in 1:nrow(params)) {
            
            eee_set <- read.delim(paste0(raw_res_path,"eee_run",params$eee_run[i],"_discarded_otus.tsv"),header=FALSE)
            eee_set <- eee_set$V1[eee_set$V1 %in% ord_tax$ASV]

            ord_counts_subset      <- ord_counts[!(rownames(ord_counts) %in% eee_set),]

            x <- abundance_filter(ord_counts_subset,
                                  cutoff = params$cutoff[i],
                                  cutoff_type = params$cutoff_type[i]
                                 )

            discarded_otus <- c(eee_set, x$discarded_clusters)

            y <- eval_res(discarded_otus, ord_tax, finbol_tax, cluster_tax, se_family, debug=FALSE)

            res <- rbind(res,
                         c(taxon=ord,
                           as.list(params[i,]),
                           as.list(y)
                          )
                        )
        }

        # Write results up to this order just in case
        write.table(res,res_file, sep="\t", row.names=FALSE)
    }

    res <- add_hexapoda_res(res, params)

    # Write final version of result table
    write.table(res,res_file, sep="\t", row.names=FALSE)
}


