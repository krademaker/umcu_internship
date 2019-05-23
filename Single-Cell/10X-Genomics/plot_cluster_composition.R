# TITLE:    plot_cluster_composition.R
# ABOUT:    Script to visualize cluster composition for 10x Genomics mouse scRNA-seq datasets (full and subsets)
# INPUT:    cluster_composition.csv: File detailing the percentage of cells (column Percentage), cluster ID (column Cluster) and dataset (column Dataset)
# OUTPUT:   cluster_composition_10x_genomics.pdf: Output plot stored in PDF file.
# AUTHOR:   Koen Rademaker
# DATE:     23 May 2019

########## Load required libraries ##########
library(ggplot2)

########## Load input data ##########
input_data <- read.csv('~/Desktop/10x_cluster_compositions.csv')

########## Plot cluster composition and save to file ##########
pdf('cluster_composition_10x_genomics.pdf', paper='a4r', width=11.69, height=8.27)
p <- ggplot(input_data, aes(Cluster, Percentage)) +
    geom_point(aes(colour = Dataset)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    scale_x_continuous('Cluster', labels = as.character(input_data$Cluster), breaks = input_data$Cluster) +
    ggtitle('Cluster composition for 10x Genomics subsets') +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab('Percentage of cells')
p
dev.off()