# Make main plot illustrating NEEAT

# Required packages: 
library(gridExtra)
library(ggplot2)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(scales)
library(ggnewscale)
library(patchwork)

source("../code/pareto_front_fxn.R")


################## Plot A #######################

# Define the plot function for the benchmark plot
plot_benchmark <- function(res, taxon) {

    res <- res[res$taxon==taxon,]
    res <- res[res$false_pos_comb<0.5 & res$false_neg_comb<1.0,]

    ggplot() +
        geom_point(data = res, size = 4, aes(x = false_neg_comb, y = false_pos_comb, shape=Method, color=Method)) +
        scale_shape_manual(values=c(16,15,3,17)) +
        scale_color_manual(values=c("blue","orange","black","green")) +
        theme_minimal() +
        theme(legend.position = "bottom") +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "inches"),
                axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0),size=16),
                axis.title.y=element_text(margin=margin(t=0, r=10, b=0, l=0),size=16),
                axis.text=element_text(size=16),
                legend.text=element_text(size=14),
                legend.title=element_text(size=16),
                plot.tag = element_text(size=22)
              ) +
        xlim(0.0, 1.0) +
        ylim(0.0, 0.5) +
        guides(shape = guide_legend()) +
        ylab("False positive rate (combined estimate)") +
        xlab("False negative rate (combined estimate)")
}

# Read in abundance results
ab_res <- read.delim("../results/abundance_res.tsv")
ab_res <- ab_res[ab_res$count_type=="raw" & ab_res$cutoff_type=="max",]
ab_res$Method <- "Abundance"

# Read in LULU results
lulu_res <- read.delim("../results/lulu_res.tsv")
lulu_res <- lulu_res[lulu_res$aligner=="vsearch",]
lulu_res$Method <- "LULU"

# Read in HMM results
hmm_res <- read.delim("../results/hmm_res.tsv")
hmm_res$Method <- "HMM"

# Read in NEEAT results
neeat_res <- read.delim("../results/eeea_res.tsv")
neeat_res <- neeat_res[neeat_res$taxon=="Hexapoda",]
neeat_res <- pareto_front(neeat_res)
neeat_res$Method <- "NEEAT"

# Merge results
cols <- c("Method","taxon","false_pos_comb","false_neg_comb")
res <- rbind(ab_res[,cols],lulu_res[,cols],hmm_res[,cols],neeat_res[,cols])

# Generate the plot
plot_A <- plot_benchmark(res,"Hexapoda")


################## Plot B #######################

# Define the plot function for the annotation plot
plot_tax_annot <- function(df, ylab, reads=FALSE) {

    if (!reads)
        df$proportion <- df$prop_count
    else
        df$proportion <- df$prop_reads

    ggplot() +
        geom_line (data=df, aes(x=rank_level, y=proportion, group=group, color=method, linetype=country)) +
        geom_point(data=df, aes(x=rank_level, y=proportion, group=group, color=method, shape=country), size=4) +
        theme_minimal() +
        theme(legend.position = "bottom") +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "inches"),
                axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0),size=16),
                axis.title.y=element_text(margin=margin(t=0, r=10, b=0, l=0),size=16),
                axis.text=element_text(size=16),
                legend.text=element_text(size=14),
                legend.title=element_text(size=16),
                plot.tag = element_text(size=22)
              ) +
        scale_x_continuous(breaks = 1:6, labels = ranks) +
        scale_color_manual(values = c("blue","orange","green")) +
        guides(
            color = guide_legend(title="Method"),   # Linetype and color legend
            shape = guide_legend(title="Country"),  # Shape legend
            linetype = guide_legend(title="Country")    # Linetype legend
        ) +
        ylab(ylab) +
        xlab("Rank")
}

# Read in and augment result data frame
annotation_res <- read.delim("../results/annotation_res.tsv")
annotation_res$group <- factor(paste(annotation_res$method,annotation_res$country,sep="."))
annotation_res$country <- factor(annotation_res$country,levels=c("Sweden","Madagascar"))
annotation_res$method <- factor(annotation_res$method)
ranks <- c("Phylum","Class","Order","Family","Genus","Species")
annotation_res$rank_level <- match(annotation_res$rank,ranks)

# Get subsets for ASVs and clusters
asv_res <- annotation_res %>% filter(count_type=="ASVs")
cluster_res <- annotation_res %>% filter(count_type=="clusters")

# Generate plot
plot_B <- plot_tax_annot(asv_res,"Proportion of ASV reads", reads=TRUE)


########### Data for plots C & D ################

# TODO: Set path to processed data
data_path <- "~/dev/figshare-repos/iba/processed_data/v2/"

# Read in data
ref_data <- read.delim("../evaluation_data/se_fauna_family.tsv")
T1 <- read.delim(paste0(data_path,"cluster_taxonomy_SE.tsv"))
T2 <- read.delim(paste0(data_path,"cleaned_noise_filtered_cluster_taxonomy_SE.tsv"))

# Get taxonomic annotation of ref ASV
T1 <- T1[T1$representative==1,]

# Limit to the shared set of families
T2 <- T2[T2$Family %in% ref_data$Family,]
T1 <- T1[T1$Family %in% T2$Family,]

# Get family counts
df_raw <- data.frame(table(T1$Family))
colnames(df_raw) <- c("Family","OTUs")

df_filtered <- data.frame(table(T2$Family))
colnames(df_filtered) <- c("Family","OTUs")

# Add data from SE fauna
df_raw <- merge(df_raw, ref_data)
df_filtered <- merge(df_filtered,ref_data)


########### Fxns for plots C & D ################

# Dummy data for plot
dummy_data <- data.frame(
    x = c(1, 10),
    y1 = c(1, 10),
    y2 = c(0.5, 5)
)

# Match plot
plot_match <- function(df, lbl_pos) {

    ggplot(df, aes(x=Known_2017, y=OTUs)) + 
        geom_point(alpha = 1.0, size = 3, aes(colour=Known_2017/Estimated)) +
        geom_text(data = subset(df, log(OTUs) > log(Known_2017) + 1.5),
                  aes(label = Family), color="black", alpha = 1.0, hjust = 0.3, vjust = -0.7) +
#       expand_limits(x = 0, y = 0) +
        scale_x_continuous(trans="log10",breaks=c(1,10,100,1000,10000,100000,1000000),
                           labels=c("1","10","100","1000","10000","100000","1000000")) +  
        scale_y_continuous(trans="log10",breaks=c(1,10,100,1000,10000,100000,1000000),
                           labels=c("1","10","100","1000","10000","100000","1000000")) +  
        scale_colour_viridis_c(option = "rocket", direction = -1, name = "Est. prop. known") +
        theme_minimal() +
        theme(legend.position = lbl_pos) +
        geom_abline(intercept = 0, slope = 1, linetype=1) +
        geom_abline(intercept = log10(0.5), slope = 1, linetype=3) +
        geom_line(data = dummy_data, aes(x = x, y = y1, linetype = "100% (y = x)")) +
        geom_line(data = dummy_data, aes(x = x, y = y2, linetype = "50% (y = 0.5x)")) +
        scale_linetype_manual(name = "Lines", values = c("100% (y = x)" = 1, "50% (y = 0.5x)" = 3)) +
        guides(linetype = guide_legend(override.aes = list(colour = c("black", "black")))) +
        coord_cartesian(clip = "off") +
        ylab("No. of clusters found") +
        xlab("No. of known species") +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 
                           "inches"), 
            axis.title.x = element_text(vjust=-1.0),
            axis.title.y = element_text(vjust=5.0),
            axis.text=element_text(size=16),
            axis.title=element_text(size=16),
            legend.text=element_text(size=12),
            legend.title=element_text(size=14),
            plot.tag = element_text(size=22)
        )
}


########### Make plots C & D ################

plot_C <- plot_match(df_raw,"none")
plot_D <- plot_match(df_filtered,"right")


########### Make final figure ###############
ggsave("Fig_neeat_main.pdf",
       width=20.0,
       height=20.0,
       plot=(plot_A + plot_B) / (plot_C + plot_D) +
            plot_annotation(tag_levels="A") +
            plot_layout(heights=c(1,1.3))
)

