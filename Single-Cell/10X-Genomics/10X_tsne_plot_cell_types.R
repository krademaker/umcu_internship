# TITLE:    10X_tsne_plot_cell_types.R
# ABOUT:    Script to plot clustering/cell types for the 10x Genomics dataset (full/108K cells/108K cell with QC applied)
# INPUT:    tsne_clusters_10x_genomics_graphclust.tsv.gz: t-SNE coordinates and cluster identities for the full dataset
# INPUT:    tsne_clusters_10x_genomics_kmeans.tsv.gz: t-SNE coordinates and cluster identities for the 108K subset
# INPUT:    tsne_cell_types_10x_genomics_108k_qc_kmeans.tsv.gz: t-SNE coordinates and cell type identities for the 108K subset with QC applied
# AUTHOR:   Koen Rademaker
# DATE:     20 June 2019

########## Load required libraries ##########
library(ggplot2)
library(dplyr)
library(ggrepel)

########## Load input data ##########
graphclust_all_cells_df <- read.delim('~/umcu_internship/Single-Cell/data/tsne_clusters_10x_genomics_graphclust.tsv.gz')
kmeans_all_cells_df <- read.delim('~/umcu_internship/Single-Cell/data/tsne_clusters_10x_genomics_kmeans.tsv.gz')
kmeans_qc_cells_df <- read.delim("~/umcu_internship/Single-Cell/data/tsne_cell_types_10x_genomics_108k_qc_kmeans.tsv.gz")

########## Format data for plotting ##########
# (1) Graph-based clustering
graphclust_all_cells_df$Cluster <- factor(graphclust_all_cells_df$Cluster)
graphclust_all_cells_labels <- data.frame(Cluster = levels(graphclust_all_cells_df$Cluster),
                                          label = levels(graphclust_all_cells_df$Cluster))
graphclust_all_cells_labels_2 <- graphclust_all_cells_df %>%
    group_by(Cluster) %>%
    summarize(TSNE.1 = mean(TSNE.1), TSNE.2 = mean(TSNE.2)) %>%
    left_join(graphclust_all_cells_labels)

# (2) K-means clustering
kmeans_all_cells_df$Cluster <- factor(kmeans_all_cells_df$Cluster)
kmeans_all_cells_labels <- data.frame(Cluster = levels(kmeans_all_cells_df$Cluster),
                                      label = levels(kmeans_all_cells_df$Cluster))
kmeans_all_cells_labels_2 <- kmeans_all_cells_df %>%
    group_by(Cluster) %>%
    summarize(TSNE.1 = mean(TSNE.1), TSNE.2 = mean(TSNE.2)) %>%
    left_join(kmeans_all_cells_labels)

# (3) K-means clustering, assigned cell types
kmeans_qc_cells_df$Cluster <- factor(kmeans_qc_cells_df$Cluster)
kmeans_qc_cells_labels <- data.frame(Cluster = levels(kmeans_qc_cells_df$Cluster), 
                                     label = levels(kmeans_qc_cells_df$Cluster))
kmeans_qc_cells_labels_2 <- kmeans_qc_cells_df %>%
    group_by(Cluster) %>%
    summarize(TSNE.1 = mean(TSNE.1), TSNE.2 = mean(TSNE.2)) %>%
    left_join(kmeans_qc_cells_labels)

########## Plot for all cells using graph-based clustering ##########
ggplot(graphclust_all_cells_df, aes(x = TSNE.1, y = TSNE.2, color = Cluster)) +
    geom_point(size = 1, alpha = 1) +
    theme_classic() +
    ggtitle('') +
    xlab('t-SNE dimension 1') +
    ylab('t-SNE dimension 2') +
    theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size=2)) +
    theme(legend.position = 'none') +
    ggrepel::geom_label_repel(data = graphclust_all_cells_labels_2, aes(label = label), size=5)

########## Plot for all cells using k-means clustering ##########
ggplot(kmeans_all_cells_df, aes(x = TSNE.1, y = TSNE.2, color = Cluster)) +
    geom_point(size = 1, alpha = 1) +
    theme_classic() +
    ggtitle('') +
    xlab('t-SNE dimension 1') +
    ylab('t-SNE dimension 2') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = 'none') +
    ggrepel::geom_label_repel(data = kmeans_all_cells_labels_2, aes(label = label), size=5)

########## Plot for 108K subset (QC applied) ##########
ggplot(kmeans_qc_cells_df, aes(x = TSNE.1, y = TSNE.2, color = Cluster)) +
    geom_point(size = 1, alpha = 1) +
    theme_classic() +
    ggtitle('') +
    xlab('t-SNE dimension 1') +
    ylab('t-SNE dimension 2') +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = 'none') +
    ggrepel::geom_label_repel(data = kmeans_qc_cells_labels_2, aes(label = label), size=4)
