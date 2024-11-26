# Evaluate the HMM results for different
# settings of the bitscore threshold
source("../code/evaluate_fxns.R")
source("../code/get_data_fxns.R")

data_path <- "../data/"

D <- read.delim("../raw_hmm_results/bitscore_filtered.txt")
D <- rbind(D,read.delim("../raw_hmm_results/bitscore_outliers.txt"))

# TODO: Set correct paths to IBA data files
raw_data_path <- "~/dev/figshare-repos/iba/raw_data/v4/"
processed_data_path <- "~/dev/figshare-repos/iba/processed_data/v2/"

# Set file names
meta_f <- paste0(raw_data_path,"CO1_sequencing_metadata_SE.tsv")
counts_f <- paste0(processed_data_path,"cluster_counts_SE.tsv")
taxonomy_f <- paste0(processed_data_path,"cluster_taxonomy_SE.tsv")

T <- get_se_cluster_taxonomy_samples_hexapoda(meta_f, counts_f, taxonomy_f)
orders <- unique(T$Order[!grepl("unclassified",T$Order) & !grepl("_X",T$Order)])

finbol_tax <- read.delim("../evaluation_data/finbol_taxonomy.tsv")
se_family <- read.delim("../evaluation_data/se_fauna_family.tsv")

res_file <- "../results/hmm_res.tsv"

bitscores <- seq(from=160,to=300,by=5)

params <- data.frame()
for (i in 1:length(bitscores))
    params <- rbind(params,list(run=i,method="hmm",cutoff=bitscores[i]))

res <- data.frame()
for (ord in orders) {

    if (sum(T$Order==ord & T$representative==1) == 1)
        next

    ord_tax <- T[T$Order==ord & T$representative==1,]

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
        write.table(res,res_file, sep="\t", row.names=FALSE)
    }
}

res <- add_hexapoda_res(res,params)

# Write final version of result table
write.table(res,res_file, sep="\t", row.names=FALSE)

