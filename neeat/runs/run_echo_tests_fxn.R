source("../code/echo_filter.R")
source("../code/evaluate_fxns.R")

run_echo_tests <- function(params) {

    data_path <- "../data/"
    raw_res_path <- "../raw_results/"

    hex_tax <- read.delim(paste0(data_path,"Hexapoda_taxonomy.tsv"), sep="\t", header=TRUE)
    orders <- unique(hex_tax$Order[!grepl("_X",hex_tax$Order) & !grepl("unclassified",hex_tax$Order)])

    res_file <- paste0("../results/echo_res_",params$run[1],".tsv")

    res <- data.frame()

    for (ord in orders) {

        cat("Processing",ord,"\n")

        ord_tax <- read.delim(paste0(data_path,ord,"_taxonomy.tsv"), sep="\t", header=TRUE)
        ord_counts <- read.delim(paste0(data_path,ord,"_cal_counts.tsv"), sep="\t", header=TRUE, row.names=1)
        ord_matchlist <- read.delim(paste0(data_path,ord,"_matchlist.tsv"), sep="\t", header=FALSE)

        # Cycle over runs for this order
        for (i in 1:nrow(params)) {

            x <- echo_filter(ord_counts,
                             ord_matchlist,
                             min_match = params$min_match[i],
                             min_overlap = params$min_overlap[i],
                             read_ratio_type = params$read_ratio_type[i],
                             max_read_ratio = params$max_read_ratio[i],
                             require_corr = params$require_corr[i])
       
            # Save essential results to files for
            # optimizing the neeat algorithm
            outfile <- paste0(raw_res_path,"echo_run",params$run[i], "_discarded_otus.tsv")
            cat(x$discarded_clusters, sep='\n', file=outfile, append=TRUE)

            y <- eval_res(x$discarded_clusters, ord_tax, debug=FALSE)

            res <- rbind(res, c(taxon=ord,
                                as.list(params[i,]),
                                as.list(y)
                                )
                        )
        }

        # Write results up to this order just in case
        write.table(res,res_file, sep="\t", row.names=FALSE)
    }

    # Compute summary result for all orders

    # Cycle over runs
    for (i in 1:nrow(params)) {

        x <- res$run==params$run[i]

        res <- rbind(res,
                     c(taxon="Hexapoda",
                       as.list(params[i,]),
                       clusters         = sum(res[x,"clusters"]),
                       removed_clusters = sum(res[x,"removed_clusters"]),
                       unique_bins      = sum(res[x,"unique_bins"]),
                       removed_uniques  = sum(res[x,"removed_uniques"]),
                       duplicate_bins   = sum(res[x,"duplicate_bins"]),
                       remaining_dups   = sum(res[x,"remaining_dups"]),
                       false_pos        = sum(res[x,"removed_uniques"])/sum(res[x,"unique_bins"]),
                       false_neg        = sum(res[x,"remaining_dups"])/sum(res[x,"duplicate_bins"])
                       )
                    )
    }

    # Write final version of result table
    write.table(res,res_file, sep="\t", row.names=FALSE)
}

