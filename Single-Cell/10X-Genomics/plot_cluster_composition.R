# TITLE:    plot_cluster_composition.R
# ABOUT:    Script to visualize cluster composition for 10x Genomics mouse scRNA-seq datasets (full and subsets)
# INPUT:    cluster_composition.csv: File detailing the percentage of cells (column Percentage), cluster ID (column Cluster) and dataset (column Dataset)
# AUTHOR:   Koen Rademaker
# DATE:     23 May 2019

########## Load required libraries ##########
library(ggplot2)

########## Load input data ##########
data <- read.csv('~/umcu_internship/Single-Cell/data/10x_Genomics_cluster_compositions_kmeans_20.csv.gz')
data$Cluster <- as.factor(data$Cluster)

########## Plot cluster composition and save to file ##########
ggplot(data, aes(Cluster, Percentage)) +
    geom_bar(aes(fill = Dataset), position =  'dodge', stat = 'identity') +
    theme_classic() +
    scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
    theme(legend.position = 'bottom', legend.title = element_text(size=12.5), legend.text = element_text(size = 10), axis.text = element_text(size = 12.5), axis.title = element_text(size = 12.5)) +
    xlab('Cluster') +
    ylab('Percentage of cells')