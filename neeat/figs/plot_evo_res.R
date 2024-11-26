# Plot summary of evo tests

# Required packages
library(gridExtra)
library(ggplot2)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(ggnewscale)
library(patchwork)

# Define a plot function for each panel
plot_res <- function(df, taxon) {

    # Get results for taxon
    df <- df[df$taxon==taxon,]
    df <- df[df$false_pos_comb<0.5,]

    ggplot() +
        geom_line (data=df, aes(x=false_neg_comb, y=false_pos_comb, group=method, color=method,linetype=method)) +
        geom_point (data=df, size=2, aes(x=false_neg_comb, y=false_pos_comb, group=method, color=method, shape=method)) +
        scale_color_manual(values=c("blue","orange","green","red","purple")) +
        theme_minimal() +
        theme(plot.margin = unit(c(0.5, 0.5, 1.5, 0.5), "cm"),
              axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0),size=14),
              axis.title.y=element_text(margin=margin(t=0, r=10, b=0, l=0),size=14)) +
        xlim(0.0, 1.0) +
        ylim(0.0, 0.5) +
        guides(
            linetype = guide_legend(title="Method"),
            color = guide_legend(title="Method"),
            shape = guide_legend(title="Method")
        ) +
        ylab("False positive rate (combined estimate)") +
        xlab("False negative rate (combined estimate)") +
        labs(title = taxon)
}

# Set data path of result files
data_path <- "../results/"

# Read in evo and hmm results
res <- read.delim(paste0(data_path, "evo_res.tsv"))
x <- character(nrow(res))
for (i in 1:nrow(res)) {
    if (res$dist_type[i]=="dadn")
        x[i] <- "dAdN."
    else
        x[i] <- "dWAdN."
    if (res$require_overlap[i])
        x[i] <- paste0(x[i],"local")
    else
        x[i] <- paste0(x[i],"global")
}
res$method <- x
hmm_res <- read.delim(paste0(data_path, "hmm_res.tsv"))
hmm_res$method <- "HMM"
cols <- c("taxon","method","false_pos_comb","false_neg_comb")
df <- rbind(res[,cols],hmm_res[,cols])

# Define taxa of interest
taxa <- c("Lepidoptera","Coleoptera","Hymenoptera","Diptera","Hemiptera","Hexapoda")

# Make each panel plot (without legend)
p1 <- plot_res(df, taxa[1])
p2 <- plot_res(df, taxa[2])
p3 <- plot_res(df, taxa[3])
p4 <- plot_res(df, taxa[4])
p5 <- plot_res(df, taxa[5])
p6 <- plot_res(df, taxa[6])

# Combine plots using patchwork with a common legend at the bottom
ggsave("Fig_evo_res.jpg", width=11.7, height=9.8,
       plot=p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(guides="collect", axes="collect") & theme(legend.position="bottom"))

