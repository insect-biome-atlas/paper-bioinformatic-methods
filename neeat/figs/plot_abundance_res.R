# Plot summary of abundance tests

# Required packages: 
library(gridExtra)
library(ggplot2)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(ggnewscale)
library(patchwork)

# Define a plot function for each panel
# with different version for different
# evaluation metrics
plot_res <- function(df, taxon) {

    df <- df[res$taxon==taxon,]
    df <- df[df$false_neg_comb<1.0 & df$false_pos_comb<0.5,]

    ggplot() +
        geom_line (data=df, aes(x=false_neg_comb, y=false_pos_comb, group=group, linetype=factor(count_type), color=factor(cutoff_type))) +
        geom_point(data=df, aes(x=false_neg_comb, y=false_pos_comb, group=group, shape=factor(cutoff_type), color=factor(cutoff_type))) +
        scale_color_manual(values=c("green","blue")) +
        theme_minimal() +
        theme(plot.margin = unit(c(0.5, 0.5, 1.5, 0.5), "cm"),
              axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0),size=14),
              axis.title.y=element_text(margin=margin(t=0, r=10, b=0, l=0),size=14)) +
        xlim(0.0, 1.00) +
        ylim(0.0, 0.5) +
        guides(
            linetype = guide_legend(title="Count type"),  # Line type legend
            shape = guide_legend(title="Cutoff type"),  # Shape legend
            color = guide_legend(title="Cutoff type")   # Color legend
        ) +
        ylab("False positive rate (combined estimate)") +
        xlab("False negative rate (combined estimate)") +
        labs(title = taxon)
}


# Set data path of result files
data_path <- "../results/"

# Read in abundance results
res <- read.delim(paste0(data_path, "abundance_res.tsv"))
res$group <- factor(paste(res$count_type,res$cutoff_type,sep="."))
res$threshold <- res$cutoff
for (i in 1:nrow(res)) { if (res$count_type[i] %in% c("sample_proportional","tot_proportional")) res$threshold[i] <- round(res$cutoff[i]*1E6) }

# Order x values correctly if they are the same
x_diff <- (1E-6 * match(res$threshold,sort(unique(res$threshold),decreasing=TRUE)))
res$false_neg <- res$false_neg + x_diff
res$false_neg_spurious <- res$false_neg_spurious + x_diff
res$false_neg_comb <- res$false_neg_comb + x_diff


# Make figure
# -----------

# Define taxa of interest
taxa <- c("Lepidoptera","Coleoptera","Hymenoptera","Diptera","Hemiptera","Hexapoda")

file_name <- paste0("Fig_abundance_res.jpg")

# Make each panel plot
p1 <- plot_res(res, taxa[1])
p2 <- plot_res(res, taxa[2])
p3 <- plot_res(res, taxa[3])
p4 <- plot_res(res, taxa[4])
p5 <- plot_res(res, taxa[5])
p6 <- plot_res(res, taxa[6])

# Combine plots using patchwork with a common legend at the bottom
ggsave(file_name, width=11.7, height=9.8,
       plot=p1 + p2 + p3 + p4 + p5 + p6 +
            plot_layout(ncol=3, guides="collect", axes="collect") &
            theme(legend.position="bottom")
      )

