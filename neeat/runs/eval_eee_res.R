# Evaluate eee results
source("../code/evaluate_fxns.R")
source("../code/get_data_fxns.R")

data_path <- "../data/"
raw_res_path <- "../raw_results/"

params <- read.delim("eee_run_params.tsv")

T <- get_se_cluster_taxonomy_samples_hexapoda()
orders <- unique(T$Order[!grepl("unclassified",T$Order) & !grepl("_X",T$Order)])

finbol_tax <- read.delim("../evaluation_data/finbol_taxonomy.tsv")
se_family <- read.delim("../evaluation_data/se_fauna_family.tsv")

res_file <- "../results/eee_res.tsv"
res <- data.frame()

for (ord in orders) {

    ord_tax <- read.delim(paste0(data_path,ord,"_taxonomy.tsv"))

    for (i in 1:nrow(params)) {

        file <- paste0(raw_res_path,"eee_run",params$run[i],"_discarded_otus.tsv")
        discarded_otus <- read.delim(file,header=FALSE)
        asv_tax <- T[T$Order==ord,]

        y <- eval_res(discarded_otus$V1, ord_tax, finbol_tax, asv_tax, se_family, debug=FALSE)

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

