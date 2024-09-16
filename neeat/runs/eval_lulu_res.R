# Evaluate the LULU results for default setting
source("../code/evaluate_fxns.R")
source("../code/get_data_fxns.R")

data_path <- "../data/"

D <- read.delim("../raw_lulu_results/lulu_res.tsv")

T <- get_se_cluster_taxonomy_samples_hexapoda()
orders <- unique(T$Order[!grepl("unclassified",T$Order) & !grepl("_X",T$Order)])

T1 <- T[T$representative==1,]
D$ASV <- T1$ASV[match(D$cluster,T1$cluster)]

finbol_tax <- read.delim("../evaluation_data/finbol_taxonomy.tsv")
se_family <- read.delim("../evaluation_data/se_fauna_family.tsv")

res_file <- "../results/lulu_res.tsv"

params <- data.frame(list(run=1,method="LULU",min_match=84,min_overlap=0.95,read_ratio_type="max",max_read_ratio=1.0))

res <- data.frame()
for (ord in orders) {

    ord_tax <- read.delim(paste0(data_path,ord,"_taxonomy.tsv"))

    discarded_otus <- D$ASV[D$numt==TRUE & D$ASV %in% ord_tax$ASV]

    cluster_tax <- T[T$Order==ord,]

    y <- eval_res(discarded_otus, ord_tax, finbol_tax, cluster_tax, se_family, debug=FALSE)

    res <- rbind(res,
                 c(taxon=ord,
                   as.list(params[1,]),
                   as.list(y)
                  )
                )

    # Write results up to this order just in case
    write.table(res,res_file, sep="\t", row.names=FALSE)
}

res <- add_hexapoda_res(res,params)

# Write final version of result table
write.table(res,res_file, sep="\t", row.names=FALSE)

