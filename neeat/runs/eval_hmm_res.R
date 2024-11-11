# Evaluate the HMM results for different
# settings of the bitscore threshold
source("../code/evaluate_fxns.R")
source("../code/get_data_fxns.R")

data_path <- "../data/"

D <- read.delim("../raw_hmm_results/bitscore_filtered.txt")
D <- rbind(D,read.delim("../raw_hmm_results/bitscore_outliers.txt"))

T <- get_se_cluster_taxonomy_samples_hexapoda()
orders <- unique(T$Order[!grepl("unclassified",T$Order) & !grepl("_X",T$Order)])

finbol_tax <- read.delim("../evaluation_data/finbol_taxonomy.tsv")
se_family <- read.delim("../evaluation_data/se_fauna_family.tsv")

res_file <- "../results/hmm_res.tsv"

bitscores <- seq(from=163,to=243,by=10)

params <- data.frame()
for (i in 1:length(bitscores))
    params <- rbind(params,list(run=i,method="hmm",cutoff=bitscores[i]))

res <- data.frame()
for (ord in orders) {
    taxfile <- paste0(data_path,ord,"_taxonomy.tsv")
    if (!file.exists(taxfile)) next
    ord_tax <- read.delim(taxfile)

    for (i in 1:nrow(params)) {

        cutoff <- params[i,"cutoff"]
    
        discarded_otus <- D$ASV[D$bitscore <= cutoff & D$ASV %in% ord_tax$ASV]
        cluster_tax <- T[T$Order==ord,]

        y <- eval_res(discarded_otus, ord_tax, finbol_tax, cluster_tax, se_family, debug=FALSE)

        res <- rbind(res,
                     c(taxon=ord,
                       as.list(params[i,]),
                       as.list(y)
                      )
                    )

        # Write results up to this order just in case
        #write.table(res,res_file, sep="\t", row.names=FALSE)
    }
}

res <- add_hexapoda_res(res,params)

# Write final version of result table
write.table(res,res_file, sep="\t", row.names=FALSE)

