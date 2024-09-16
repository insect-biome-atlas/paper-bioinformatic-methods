source("../code/evo_filter.R")
source("../code/evaluate_fxns.R")

run_eee_tests <- function(params) {

    data_path <- "../data/"
    raw_res_path <- "../raw_results/"

    hex_tax <- read.delim(paste0(data_path,"Hexapoda_taxonomy.tsv"), sep="\t", header=TRUE)
    orders <- unique(hex_tax$Order[!grepl("_X",hex_tax$Order) & !grepl("unclassified",hex_tax$Order)])

    res_file <- paste0("../results/eee_res_",params$run[1],".tsv")
    res <- data.frame()

    for (ord in orders) {

        cat("Processing",ord,"\n")

        ord_tax <- read.delim(paste0(data_path,ord,"_taxonomy.tsv"), sep="\t", header=TRUE)

        ord_counts <- read.delim(paste0(data_path,ord,"_cal_counts.tsv"), sep="\t", header=TRUE, row.names=1)
        ord_evodistlist <- read.delim(paste0(data_path,ord,"_evodistlist.tsv"))

        # Cycle over runs for this order
        for (i in 1:nrow(params)) {
            
            ee_set <- read.delim(paste0(raw_res_path,"ee_run",params$ee_run[i],"_discarded_otus.tsv"),header=FALSE)
            ee_set <- ee_set$V1[ee_set$V1 %in% ord_tax$ASV]

            ord_counts_subset      <- ord_counts[!(rownames(ord_counts) %in% ee_set),]
            ord_evodistlist_subset <- ord_evodistlist[!(ord_evodistlist$asv1 %in% ee_set) & !(ord_evodistlist$asv2 %in% ee_set),]

            x <- evo_filter(ord_counts_subset,
                            ord_evodistlist_subset,
                            require_overlap = params$require_overlap[i],
                            dist_type       = params$dist_type[i],
                            dist_threshold  = params$dist_threshold[i]
                           )

            discarded_otus <- c(ee_set, x$discarded_clusters)

            # Save essential results to files for use
            # in optimizing the neeat algorithm
            outfile <- paste0(raw_res_path,"eee_run",params$run[i], "_discarded_otus.tsv")
            cat(discarded_otus, sep='\n', file=outfile, append=TRUE)

            y <- eval_res(discarded_otus, ord_tax, debug=FALSE)

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

