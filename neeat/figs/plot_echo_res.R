# Plot summary of echo tests

# Required packages: 
library(gridExtra)
library(ggplot2)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(ggnewscale)
library(patchwork)

# Define plot function for each panel
plot_res <- function(res, taxon) {

    # Get results for taxon
    res <- res[res$taxon==taxon,]

    ggplot() +
        geom_line (data=res, aes(x=false_neg_comb, y=false_pos_comb, group=Case, color=Case,
                                 linetype=Case)) +
        geom_point(data=res, size=2, aes(x=false_neg_comb, y=false_pos_comb, group=Case, color=Case,
                                 shape=Case)) +
        scale_shape_manual(values=c(16,15,17,18,3)) +
        scale_color_manual(name="Method", values=c("Red","Orange","Green","Blue","Black")) +
        theme_minimal() +
        theme(plot.margin=unit(c(0.5, 0.5, 1.5, 0.5), "cm"),
              axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0),size=14),
              axis.title.y=element_text(margin=margin(t=0, r=10, b=0, l=0),size=14)) +
        xlim(0.0, 1.0) +
        ylim(0.0, 0.10) +
        guides(
            shape = guide_legend(title="Method"),
            linetype = guide_legend(title="Method"),
            color = guide_legend(title="Method")
        ) +
        ylab("False positive rate (combined estimate)") +
        xlab("False negative rate (combined estimate)") +
        labs(title = taxon)
}



# Set data path of result files
data_path <- "../results/"

# Read in echo and lulu results
res <- read.delim(paste0(data_path, "echo_res.tsv"))
lulu_res <- read.delim(paste0(data_path, "lulu_res.tsv"))
lulu_res <- lulu_res[lulu_res$aligner=="vsearch",]

# Add case labels
res$Case <- ""
res$Case[res$read_ratio_type=="max" & res$require_corr==TRUE] <- "Corr.max"
res$Case[res$read_ratio_type=="mean" & res$require_corr==TRUE] <- "Corr.mean"
res$Case[res$read_ratio_type=="max" & res$require_corr==FALSE] <- "No.corr.max"
res$Case[res$read_ratio_type=="mean" & res$require_corr==FALSE] <- "No.corr.mean"

lulu_res$Case <- "LULU"
lulu_res$max_read_ratio <- 1.0
lulu_res$require_corr <- FALSE

# Merge
colnames <- c("taxon","Case","max_read_ratio","read_ratio_type","require_corr","false_pos","false_neg","false_pos_finbol","false_neg_spurious","false_pos_comb","false_neg_comb")
res <- rbind(res[,colnames], lulu_res[,colnames])
res$Case <- factor(res$Case, levels = c("Corr.max","Corr.mean", "No.corr.max", "No.corr.mean", "LULU"))


# Generate plot
# ----------------

# Define taxa of interest
taxa <- c("Lepidoptera","Coleoptera","Hymenoptera","Diptera","Hemiptera","Hexapoda")

# Make figure
file_name <- paste0("Fig_echo_res.jpg")

# Make each panel plot
p1 <- plot_res(res, taxa[1])
p2 <- plot_res(res, taxa[2])
p3 <- plot_res(res, taxa[3])
p4 <- plot_res(res, taxa[4])
p5 <- plot_res(res, taxa[5])
p6 <- plot_res(res, taxa[6])

# Combine plots using patchwork with a common legend at the bottom
ggsave(file_name, width=11.7, height=8.5,
       plot=p1 + p2 + p3 + p4 + p5 + p6 +
            plot_layout(ncol=3, guides="collect", axes="collect") &
            theme(legend.position="bottom")
      )


