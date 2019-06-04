# TITLE:    10X_tsne_plot_cell_types.R
# ABOUT:    Script to plot clustering/cell types for the full 10x Genomics dataset and for the 108K cells subset with QC applied
# INPUT:    tsne_clusters_10x_genomics.xlsx: t-SNE coordinates and clustering results for the full dataset.
# INPUT:    tsne_cell_types_10x_genomics_108k_subset_qc.xlsx: t-SNE coordinates and cell types identities for 108K subset with QC applied
# AUTHOR:   Koen Rademaker
# DATE:     4 June 2019

########## Load required libraries ##########
library(ggplot2)
library(dplyr)
library(ggrepel)
library(readxl)

########## Load input data ##########
all_cells_df <- read_excel('~/umcu_internship/Single-Cell/data/tsne_clusters_10x_genomics.xlsx')
qc_cells_df <- read_excel('~/umcu_internship/Single-Cell/data/tsne_cell_types_10x_genomics_108k_subset_qc.xlsx')

########## Format data for plotting ##########
all_cells_df$Cluster <- factor(all_cells_df$Cluster)
all_cells_labels <- data.frame(Cluster = levels(all_cells_df$Cluster),
                               label = paste0('Cluster ', levels(all_cells_df$Cluster)))
all_cells_labels_2 <- all_cells_df %>%
    group_by(Cluster) %>%
    summarize(TSNE_1 = mean(TSNE_1), TSNE_2 = mean(TSNE_2)) %>%
    left_join(all_cells_labels)


qc_cells_df$Cluster <- factor(qc_cells_df$Cluster)
qc_cells_labels <- data.frame(Cluster = levels(qc_cells_df$Cluster), 
                              label = levels(qc_cells_df$Cluster))
qc_cells_labels_2 <- qc_cells_df %>%
    group_by(Cluster) %>%
    summarize(TSNE_1 = mean(TSNE_1), TSNE_2 = mean(TSNE_2)) %>%
    left_join(qc_cells_labels)

########## Plot for all cells ##########
ggplot(all_cells_df, aes(x = TSNE_1, y = TSNE_2, color = Cluster)) +
    geom_point(size = 1, alpha = 1) +
    theme_classic() +
    ggtitle('K-means clustering (k=20) for 10x Genomics dataset') +
    xlab('t-SNE dimension 1') +
    ylab('t-SNE dimension 2') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = 'none') +
    ggrepel::geom_label_repel(data = all_cells_labels_2, aes(label = label))

########## Plot for 108K subset (QC applied) ##########
ggplot(qc_cells_df, aes(x = TSNE_1, y = TSNE_2, color = Cluster)) +
    geom_point(size = 1, alpha = 1) +
    theme_classic() +
    ggtitle('Cell types for 10x Genomics 108K subset with QC applied') +
    xlab('t-SNE dimension 1') +
    ylab('t-SNE dimension 2') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = 'none') +
    ggrepel::geom_label_repel(data = qc_cells_labels_2, aes(label = label))