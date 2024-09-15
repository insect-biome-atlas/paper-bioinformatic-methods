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

    res_file <- paste0("../results/eeea_res_",params$run[1],".tsv")

    res <- data.frame(c(taxon=character(),
                        as.list(params[FALSE,]),
                        clusters=numeric(),
                        removed_clusters=numeric(),
                        unique_bins=numeric(),
                        removed_uniques=numeric(),
                        duplicate_bins=numeric(),
                        remaining_dups=numeric(),
                        false_pos=numeric(),
                        false_neg=numeric(),
                        finbol_clusters=numeric(),
                        removed_finbol_clusters=numeric(),
                        spp_in_known=numeric(),
                        found_spp_in_known=numeric(),
                        clusters_in_known=numeric(),
                        removed_clusters_in_known=numeric(),
                        false_pos_finbol=numeric(),
                        false_neg_spurious=numeric(),
                        false_pos_comb=numeric(),
                        false_neg_comb=numeric())
                     )

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

    # Compute summary result for all hexapod orders

    # Cycle over runs
    for (i in 1:nrow(params)) {

        x <- res$run==params$run[i]

        # Do the complicated math first
        false_pos = sum(res[x,"removed_uniques"])/sum(res[x,"unique_bins"])
        false_neg = sum(res[x,"remaining_dups"])/sum(res[x,"duplicate_bins"])
        false_pos_finbol = sum(res[x,"removed_finbol_clusters"])/sum(res[x,"finbol_clusters"])
        spurious = sum(res[x,"clusters_in_known"]) - sum(res[x,"found_spp_in_known"])
        remaining_spurious = spurious - sum(res[x,"removed_clusters_in_known"])
        false_neg_spurious = remaining_spurious / spurious
        false_pos_comb = (false_pos + false_pos_finbol)/2.0
        false_neg_comb = (false_neg + false_neg_spurious)/2.0
        
        res <- rbind(res,
                     c(taxon="Hexapoda",
                       as.list(params[i,]),
                       clusters                     = sum(res[x,"clusters"]),
                       removed_clusters             = sum(res[x,"removed_clusters"]),
                       unique_bins                  = sum(res[x,"unique_bins"]),
                       removed_uniques              = sum(res[x,"removed_uniques"]),
                       duplicate_bins               = sum(res[x,"duplicate_bins"]),
                       remaining_dups               = sum(res[x,"remaining_dups"]),
                       false_pos                    = false_pos,
                       false_neg                    = false_neg,
                       finbol_clusters              = sum(res[x,"finbol_clusters"]),
                       removed_finbol_clusters      = sum(res[x,"removed_finbol_clusters"]),
                       spp_in_known                 = sum(res[x,"spp_in_known"]),
                       found_spp_in_known           = sum(res[x,"found_spp_in_known"]),
                       clusters_in_known            = sum(res[x,"clusters_in_known"]),
                       removed_clusters_in_known    = sum(res[x,"removed_clusters_in_known"]),
                       false_pos_finbol             = false_pos_finbol,
                       false_neg_spurious           = false_neg_spurious,
                       false_pos_comb               = false_pos_comb,
                       false_neg_comb               = false_neg_comb
                      )
                    )
    }

    # Write final version of result table
    write.table(res,res_file, sep="\t", row.names=FALSE)
}

