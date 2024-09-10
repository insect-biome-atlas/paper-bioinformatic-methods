# Merge taxonomy files
data_path <- "../../data/entognatha_outgroup/"
T1 <- read.delim(paste0(data_path,"collembola_taxa.csv"),sep=";")
T2 <- read.delim(paste0(data_path,"entognatha_taxa.csv"),sep=";")
T3 <- read.delim(paste0(data_path,"root_taxa.csv"),sep=";")

T1$Kingdom <- "Animalia"
T1$Phylum <- "Arthropoda"

T2$Kingdom <- "Animalia"
T2$Phylum <- "Arthropoda"

cols <- c("TipLabel","Kingdom","Phylum","Class","Order","Family","Genus","Species")

T <- rbind(T1[,cols],T2[,cols],T3[,cols])

write.table(T,paste0(data_path,"entognatha_outgroup_taxonomy.tsv"),sep="\t",row.names=FALSE)

