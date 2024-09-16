source("../code/evo_filter.R")
source("../code/evaluate_fxns.R")

run_evo_tests <- function(params) {

    data_path <- "../data/"
    raw_res_path <- "../raw_results/"

    hex_tax <- read.delim(paste0(data_path,"Hexapoda_taxonomy.tsv"), sep="\t", header=TRUE)
    orders <- unique(hex_tax$Order[!grepl("_X",hex_tax$Order) & !grepl("unclassified",hex_tax$Order)])

    res_file <- paste0("../results/evo_res_",params$run[1],".tsv")
    res <- data.frame()

    for (ord in orders) {

        cat("Processing",ord,"\n")

        ord_tax <- read.delim(paste0(data_path,ord,"_taxonomy.tsv"), sep="\t", header=TRUE)
        ord_counts <- read.delim(paste0(data_path,ord,"_cal_counts.tsv"), sep="\t", header=TRUE, row.names=1)
        ord_evodistlist <- read.delim(paste0(data_path,ord,"_evodistlist.tsv"))

        # Cycle over runs for this order
        for (i in 1:nrow(params)) {

            x <- evo_filter(ord_counts,
                            ord_evodistlist,
                            require_overlap = params$require_overlap[i],
                            dist_type       = params$dist_type[i],
                            dist_threshold  = params$dist_threshold[i]
                           )
       
            # Save essential results to files for possible use
            # in optimizing the neeat algorithm (just in case...)
            outfile <- paste0(raw_res_path,"evo_run",params$run[i], "_discarded_otus.tsv")
            cat(x$discarded_clusters, sep='\n', file=outfile, append=TRUE)

            y <- eval_res(x$discarded_clusters, ord_tax, debug=FALSE)

            res <- rbind(res, c(taxon=ord,
                                as.list(params[i,]),
                                clusters=y$clusters,
                                removed_clusters=y$removed_clusters,
                                unique_bins=y$unique_bins,
                                removed_uniques=y$removed_uniques,
                                duplicate_bins=y$duplicate_bins,
                                remaining_dups=y$remaining_dups,
                                false_pos=y$false_pos,
                                false_neg=y$false_neg
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

