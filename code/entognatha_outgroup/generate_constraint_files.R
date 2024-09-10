# Read in functions that generate hard constratints and
# partial constraints from an input tree
source("../constraint_fxns.R")

# Set working dir
old_wd <- getwd()
setwd("../../data/entognatha_outgroup/")

# Generate collembola constraints
gen_mb_con_file("collembola_tree.nwk", "collembola_constraints.nex")

# Generate higher constraints
out_file <- "higher_constraints.nex"

# Read in metadata
root_taxa <- read.delim("root_taxa.csv",sep=";")
collembola_taxa <- read.delim("collembola_taxa.csv",sep=";")
entognatha_taxa <- read.delim("entognatha_taxa.csv",sep=";")

# Print header to output file
cat("#NEXUS\n\nbegin mrbayes;\n",file=out_file)

# Assemble all tip labels and higher classification info
taxa <- data.frame(TipLabel=c(root_taxa$TipLabel, collembola_taxa$TipLabel, entognatha_taxa$TipLabel))
n <- nrow(taxa) - nrow(root_taxa)
taxa$Kingdom <- c(root_taxa$Kingdom, rep("Animalia",n))
taxa$Phylum <- c(root_taxa$Phylum, rep("Arthropoda",n))
taxa$Class <- c(root_taxa$Class, collembola_taxa$Class, entognatha_taxa$Class)
taxa$Order <- c(root_taxa$Order, collembola_taxa$Order, entognatha_taxa$Order)

# Output informative constraint partitions for higher taxa
for (kingdom in unique(taxa$Kingdom)) {
    ingroup <- taxa$TipLabel[taxa$Kingdom==kingdom]
    add_constraint(kingdom, ingroup, "higher_constraints.nex")
}

for (phylum in unique(taxa$Phylum)) {
    ingroup <- taxa$TipLabel[taxa$Phylum==phylum]
    add_constraint(phylum, ingroup, "higher_constraints.nex")
}


for (class in unique(taxa$Class)) {
    ingroup <- taxa$TipLabel[taxa$Class==class]
    add_constraint(class, ingroup, "higher_constraints.nex")
}

for (order in unique(taxa$Order)) {
    ingroup <- taxa$TipLabel[taxa$Order==order]
    add_constraint(order, ingroup, "higher_constraints.nex")
}

# Add special constraints
# -----------------------

# Opisthokonta
ingroup <- taxa$TipLabel[taxa$Kingdom!="Diphoda"]
add_constraint("Opisthokonta", ingroup, "higher_constraints.nex")

# Protostomia (Deuterostomia only includes Chordata, so not needed)
ingroup <- taxa$TipLabel[taxa$Kingdom=="Animalia" & !(taxa$Phylum=="Chordata")]
add_constraint("Protostomia", ingroup, "higher_constraints.nex")

# Lophotrochozoa
ingroup <- taxa$TipLabel[taxa$Phylum=="Mollusca" | taxa$Phylum=="Annelida"]
add_constraint("Lophotrochozoa", ingroup, "higher_constraints.nex")

# Ecdysozoa
ingroup <- taxa$TipLabel[taxa$Phylum=="Arthropoda" | taxa$Phylum=="Tardigrada" | taxa$Phylum=="Nematoda"]
add_constraint("Ecdysozoa", ingroup, "higher_constraints.nex")

# Panarthropoda
ingroup <- taxa$TipLabel[taxa$Phylum=="Arthropoda" | taxa$Phylum=="Tardigrada"]
add_constraint("Panarthropoda", ingroup, "higher_constraints.nex")

# Mandibulata
ingroup <- taxa$TipLabel[taxa$Phylum=="Arthropoda" & !(taxa$Class=="Arachnida")]
add_constraint("Mandibulata", ingroup, "higher_constraints.nex")

# Pancrustacea
ingroup <- taxa$TipLabel[taxa$Phylum=="Arthropoda" & !(taxa$Class=="Arachnida") & !(taxa$Class=="Diplopoda")]
add_constraint("Pancrustacea", ingroup, "higher_constraints.nex")

# Hexapoda
ingroup <- taxa$TipLabel[taxa$Class=="Insecta" | taxa$Class=="Collembola" | taxa$Class=="Diplura" | taxa$Class=="Protura"]
add_constraint("Hexapoda", ingroup, "higher_constraints.nex")

# Dicondylia
ingroup <- taxa$TipLabel[taxa$Order=="Zygentoma" | taxa$Order=="Ephemeroptera" | taxa$Order=="Diptera"]
add_constraint("Dicondylia", ingroup, "higher_constraints.nex")

# Pterygota
ingroup <- taxa$TipLabel[taxa$Order=="Ephemeroptera" | taxa$Order=="Diptera"]
add_constraint("Pterygota", ingroup, "higher_constraints.nex")

# Print tail to output file
cat ("end;\n", file=out_file, append=TRUE)

# Reset working directory
setwd(old_wd)

