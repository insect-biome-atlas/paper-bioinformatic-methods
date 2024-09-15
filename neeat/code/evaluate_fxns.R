# Evaluate results from cluster (=otu) filtering
# Note that discarded OTUs are assumed to be given as rep ASVs
# and not cluster names 
eval_res <- function(discarded_otus, otu_taxonomy, all_tax, se_family, finbol_tax, cutoff=0.95, debug=TRUE) {
  
  if (debug) {
    cat ("The taxonomy has", nrow(otu_taxonomy), "clusters in total\n")
    cat (sum(grepl(" ",otu_taxonomy$Species)), "clusters have rep asv annotated to species\n")
    cat ("There are", length(unique(otu_taxonomy$Species[grepl(" ",otu_taxonomy$Species)])), "unique species annotations\n")
    cat (sum(grepl("BOL",otu_taxonomy$BOLD_bin)), "clusters have rep asv annotated to BOLD bin\n")
    cat ("There are", length(unique(otu_taxonomy$BOLD_bin[grepl("BOL",otu_taxonomy$BOLD_bin)])), "unique BOLD bin annotations\n")
    cat ("There are", length(discarded_otus), "clusters removed\n")
  }

  unique_bins <- length(unique(otu_taxonomy$BOLD_bin[grepl("BOL", otu_taxonomy$BOLD_bin)]))
  remaining_unique_bins <- length(unique(otu_taxonomy$BOLD_bin[grepl("BOL", otu_taxonomy$BOLD_bin) & !(otu_taxonomy$ASV %in% discarded_otus)]))
  removed_uniques <- unique_bins - remaining_unique_bins
  false_pos <- removed_uniques / unique_bins
  if (debug) {
    cat ("False positive rate (unique  BOLD bins removed): ", removed_uniques,
         "of", unique_bins, "=", false_pos, "\n")
  }

  finbol_tax <- read.delim("../finbol/finbol_taxonomy.tsv")
  finbol_tax <- finbol_tax[finbol_tax$cluster %in% otu_taxonomy$cluster,]
  n_finbol_clusters <- length(unique(finbol_tax$cluster))
  discarded_finbol_clusters <- unique(finbol_tax$cluster[match(discarded_otus,finbol_tax$rep_asv)])
  n_discarded_finbol_clusters <- length(discarded_finbol_clusters)
  false_pos_finbol <- n_discarded_finbol_clusters / n_finbol_clusters
  if (debug) {
    cat ("False positive rate (FinBOL clusters removed): ", n_discarded_finbol_clusters,
         "of", n_finbol_clusters, "=", false_pos_finbol, "\n")
  }

  dups <- length(unique(otu_taxonomy$BOLD_bin[grepl("BOL", otu_taxonomy$BOLD_bin) & duplicated(otu_taxonomy$BOLD_bin)]))
  BOLD_bins <- otu_taxonomy[grepl("BOL", otu_taxonomy$BOLD_bin), ]
  retained_BOLD_bins <- BOLD_bins[!BOLD_bins$ASV %in% discarded_otus, ]
  remaining_dups <- length(unique(retained_BOLD_bins$BOLD_bin[duplicated(retained_BOLD_bins$BOLD_bin)]))
  false_neg = remaining_dups/dups
  if (debug) {
      cat ("False negative rate (duplicated BOLD bins remaining): ", remaining_dups, "of", dups, "=", false_neg, "\n")
  }

  se_family <- read.delim("../se_fauna/se_fauna_family.tsv")
  se_known_fams <- se_family$Family[se_family$prop_known > cutoff]
  T <- all_tax[all_tax$Family %in% se_known_fams,]
  known_spp <- sum(se_family$Known_2017[se_family$Family %in% se_known_fams])
  found_spp <- length(unique(T$Species[grepl(" ",T$Species)]))
  
  clusters_in_known <- unique(otu_taxonomy$cluster[otu_taxonomy$Family %in% se_known_fams])
  n_clusters_in_known <- length(clusters_in_known)
  removed_clusters <- otu_taxonomy$cluster[match(discarded_otus,otu_taxonomy$ASV)]
  n_removed_clusters_in_known <- sum(removed_clusters %in% clusters_in_known)
  n_remaining_clusters_in_known <- n_clusters_in_known - n_removed_clusters_in_known
  false_neg_spurious <- (n_remaining_clusters_in_known - found_spp) / (n_clusters_in_known - found_spp)

  if (debug) {
    cat ("Number of clusters in known: ", n_clusters_in_known, "\n")
    cat ("False negative rate (remaining fraction of superfluous clusters): ", n_remaining_clusters_in_known, "clusters for",
         found_spp, "found spp =", false_neg_spurious, "fraction remaining spurious clusters\n")
  }

  list(clusters=nrow(otu_taxonomy),
       removed_clusters=length(discarded_otus),
       unique_bins=unique_bins,
       removed_uniques=removed_uniques,
       duplicate_bins=dups,
       remaining_dups=remaining_dups,
       false_pos=false_pos,
       false_neg=false_neg,
       spp_in_known=known_spp,
       found_spp_in_known=found_spp,
       clusters_in_known=n_clusters_in_known,
       removed_clusters_in_known=n_removed_clusters_in_known,
       false_pos_finbol=false_pos_finbol,
       false_neg_spurious=false_neg_spurious)
}

