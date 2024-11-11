library(ape)
library(phytools)
source("../seq_fxns.R")

old_wd <- getwd()
setwd("../../data/chesters_tree/")

# Read in taxonomy for Chesters tree.
# The taxonomy has Kingdom=Eukaryota but this is a superkingdom in TOL,
# which is the taxonomy used by GBIF. Therefore, we correct the Kingdom here.
chesters_taxonomy <- read.delim("chesters_taxonomy.tsv")
chesters_taxonomy$Kingdom <- "Animalia"

# Read in Chesters tree (the ultrametric tree)
chesters_tree <- read.tree("species_level_tree.ultrametric.nwk")

# Drop tips that have data issues and are omitted from
# the taxonomy file for that reason.
missing_tips <- chesters_tree$tip.label[!(chesters_tree$tip.label %in% chesters_taxonomy$TipLabel)]
chesters_tree <- drop.tip(chesters_tree,missing_tips)

# Drop tips for which sequences have issues. Update sequences and taxonomy file
seqs <- read.FASTA("chesters_CO1.fasta")
seqs <- seqs[names(seqs) %in% chesters_taxonomy$TipLabel]
x1 <- no_stop_codons(seqs,1)
x2 <- no_stop_codons(seqs,2)
x3 <- no_stop_codons(seqs,3)
chesters_taxonomy$Start <- -1
chesters_taxonomy$Start[match(names(seqs[x3]),chesters_taxonomy$TipLabel)] <- 3
chesters_taxonomy$Start[match(names(seqs[x2]),chesters_taxonomy$TipLabel)] <- 2 
chesters_taxonomy$Start[match(names(seqs[x1]),chesters_taxonomy$TipLabel)] <- 1
problem_taxa <- (names(seqs))[!(x1|x2|x3)]
seqs <- seqs[x1|x2|x3]
chesters_taxonomy <- chesters_taxonomy[chesters_taxonomy$TipLabel %in% names(seqs),]
chesters_taxonomy$Stop <- -1
for (i in 1:length(seqs)) {
    idx <- match(names(seqs)[i],chesters_taxonomy$TipLabel)
    len <- length(seqs[[i]]) - chesters_taxonomy$Start[idx] + 1
    chesters_taxonomy$Stop[idx] <- length(seqs[[i]]) - (len %% 3)
}
chesters_tree <- drop.tip(chesters_tree,problem_taxa)
for (i in 1:nrow(chesters_taxonomy)) {
    if ((1 + chesters_taxonomy$Stop[i] - chesters_taxonomy$Start[i]) %% 3 != 0)
        cat("Problems with length of taxon", chesters_taxonomy$TipLabel[i],"\n")
}

# Replace Entognatha and outgroup part of tree
# --------------------------------------------

# Drop outgroup tips
chesters_outgroups <- chesters_taxonomy$TipLabel[chesters_taxonomy$Class!="Insecta"]
chesters_tree <- drop.tip(chesters_tree,chesters_outgroups)

# Read new outgroup tree (assuming a single tree in the file)
root_tree <- read.nexus("../entognatha_outgroup/entognatha_outgroup.mb.last_tree.tre")

# Get taxon info on root tree
entognatha_outgroup_taxonomy <- read.delim("../entognatha_outgroup/entognatha_outgroup_taxonomy.tsv")
insecta <- entognatha_outgroup_taxonomy$TipLabel[entognatha_outgroup_taxonomy$Class=="Insecta"]
diplura <- entognatha_outgroup_taxonomy$TipLabel[entognatha_outgroup_taxonomy$Class=="Diplura"]

insect_node <- getMRCA(root_tree, insecta)
diplura_node <- getMRCA(root_tree, diplura)
insert_node <- getMRCA(root_tree, c(insecta, diplura))
depths <- node.depth.edgelength(root_tree)
root_tree_height <- depths[1]     # Total "height" of tree
insect_height <- depths[insect_node]
diplura_height <- depths[diplura_node]
insert_height <- depths[insert_node]
# Insert position also computable as diplura_height - insert_height
# Here we get the appropriate branch length instead
insert_pos <- root_tree$edge.length[which(root_tree$edge[,2]==diplura_node)] 

# Rescale the chesters tree as donor tree
donor_tree <- chesters_tree
donor_tree_height <- node.depth.edgelength(donor_tree)[1]
scale_factor <- (root_tree_height - insect_height) / donor_tree_height
donor_tree$edge.length <- donor_tree$edge.length * scale_factor
donor_tree$root.edge <- insect_height - insert_height   # Set the length of the root branch to get an ultrametric tree

# Now we are ready to combine the trees
receptor_tree <- drop.tip(root_tree, insecta)
insert_node <- getMRCA(receptor_tree, diplura)  # Node number of diplura node might have changed
expanded_tree <- bind.tree(receptor_tree, donor_tree, where=insert_node, position=insert_pos)

# Expand the combined tree as the chesters lengths are presumably more accurate
# because of increased taxon sampling
expanded_tree$edge.length <- expanded_tree$edge.length / scale_factor

# Write combined tree
cat("Writing expanded tree\n")
write.tree(expanded_tree,"chesters_new_outgroups.nwk")

# Collect sequences matching this tree
# The Chesters sequences are not trimmed to the coding part yet, so this is done here
# The other sequences have already been processed in this way
seqs <- extract_coding(seqs,chesters_taxonomy)
seqs <- seqs[!(names(seqs) %in% c(chesters_outgroups))]
write.FASTA(seqs,"chesters_new_outgroups.fasta")
seqs <- read.FASTA("../entognatha_outgroup/entognatha_outgroup.fasta")
seqs <- seqs[!(names(seqs) %in% insecta)]
write.FASTA(seqs,"chesters_new_outgroups.fasta",append=TRUE)

# Collect taxonomy info for this tree (note we also want to replace braconid data with more detailed data)
chesters_drop_taxa <- c(chesters_outgroups)
chesters_taxonomy <- chesters_taxonomy[!(chesters_taxonomy$TipLabel %in% chesters_drop_taxa),]
entognatha_outgroup_taxonomy <- entognatha_outgroup_taxonomy[!(entognatha_outgroup_taxonomy$TipLabel %in% insecta),]

cols <- c("TipLabel","Kingdom","Phylum","Class","Order","Family","Genus","Species")
T <- rbind(chesters_taxonomy[,cols],entognatha_outgroup_taxonomy[,cols])
write.table(T,"chesters_new_outgroups_taxonomy.tsv",sep="\t",row.names=FALSE)

# Reset working dir
setwd(old_wd)

