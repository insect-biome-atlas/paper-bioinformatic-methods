library(dplyr)
library(ggplot2)
library(patchwork)

plot_annotation_res <- function(df, ylab, reads=FALSE) {

    if (!reads)
        df$proportion <- df$prop_count
    else
        df$proportion <- df$prop_reads

    ggplot() +
        geom_line (data=df, aes(x=rank_level, y=proportion, group=group, color=method, linetype=country)) +
        geom_point(data=df, aes(x=rank_level, y=proportion, group=group, color=method, shape=country), size=3) +
        theme_minimal() +
        theme(plot.margin = unit(c(0.5, 0.5, 1.5, 0.5), "cm"),
              axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0),size=12),
              axis.title.y=element_text(margin=margin(t=0, r=10, b=0, l=0),size=12)) +
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

# Generate plots
p1 <- plot_annotation_res(asv_res,"Proportion of ASVs")
p2 <- plot_annotation_res(asv_res,"Proportion of ASV reads", reads=TRUE)
p3 <- plot_annotation_res(cluster_res,"Proportion of clusters")
p4 <- plot_annotation_res(cluster_res,"Proportion of cluster reads", reads=TRUE)

ggsave("Fig_annotation_res.pdf", width=10.0, height=10.0, plot=(p1 + p2) / (p3 + p4) + plot_layout(guides="collect"))
