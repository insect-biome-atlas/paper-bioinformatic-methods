# Evaluate the LULU results for default setting
source("../code/evaluate_fxns.R")
source("../code/get_data_fxns.R")

data_path <- "../data/"

T <- get_se_cluster_taxonomy_samples_hexapoda()
orders <- unique(T$Order[!grepl("unclassified",T$Order) & !grepl("_X",T$Order)])

T1 <- T[T$representative==1,]


finbol_tax <- read.delim("../evaluation_data/finbol_taxonomy.tsv")
se_family <- read.delim("../evaluation_data/se_fauna_family.tsv")

res_file <- "../results/lulu_res.tsv"

params <- data.frame(list(run=c(1,3),method="LULU",min_match=84,minimum_relative_cooccurence=0.95, min_overlap=0.95,read_ratio_type="max",max_read_ratio=1.0, aligner=c("vsearch","blastn"), max_target_seqs=10, minimum_ratio=1, minimum_ratio_type="min"))

res <- data.frame()
i <- 1
for (run in params$run) {
  raw_file <- paste0("../raw_lulu_results/lulu_raw_",run,".tsv")
    if (!file.exists(raw_file)) {
      next
    }
  run_res <- data.frame()
  for (ord in orders) {
    D <- read.delim(raw_file)
    D$ASV <- T1$ASV[match(D$cluster,T1$cluster)]
    taxfile <- paste0(data_path,ord,"_taxonomy.tsv")
    if (!file.exists(taxfile)) {
      next
    }
    ord_tax <- read.delim(paste0(data_path,ord,"_taxonomy.tsv"))

    discarded_otus <- D$ASV[D$numt==TRUE & D$ASV %in% ord_tax$ASV]

    cluster_tax <- T[T$Order==ord,]

    y <- eval_res(discarded_otus, ord_tax, finbol_tax, cluster_tax, se_family, debug=FALSE)

    run_res <- rbind(run_res,
                c(taxon=ord,
                  as.list(params[i,]),
                  as.list(y)
                  )
                )
  }
  run_res <- add_hexapoda_res(run_res,params[i,])
  # Write results up to this order just in case
  #write.table(res,res_file, sep="\t", row.names=FALSE)
  res <- rbind(res,run_res)
  i <- i+1
}


# Write final version of result table
write.table(res,res_file, sep="\t", row.names=FALSE)

