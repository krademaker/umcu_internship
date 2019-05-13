# TITLE:    cluster_composition.R
# ABOUT:    Script to determine and visualize cluster composition for the full or subsampled 10x Genomics mouse brain cell dataset.
# INPUT:    [N]_graphcluster_clusters.csv_composition: File(s) detailing number of cells per cluster (N) and cluster ID (ID)
# AUTHOR:   Koen Rademaker
# DATE:     13 May 2019


########## All cells (n=1,306,217) ##########
clus_comp_full <- read.table('/home/koen/Data/SingleCell_Data/10X_Genomics/cluster_composition/1m_graphclust_clusters.csv_composition', quote='\'', comment.char='', col.names = c('N', 'ID'))
clus_comp_full <- clus_comp_full[!(clus_comp_full$ID=='Cluster'),]
clus_comp_full$ID <- as.numeric(clus_comp_full$ID)
clus_comp_full$percentage = clus_comp_full$N/sum(clus_comp_full$N)
clus_comp_full <- clus_comp_full[order(clus_comp_full$ID),]
clus_comp_full$percentage <- clus_comp_full$percentage*100


########## 1/2 of cells (n=653,108) ##########
clus_comp_half <- read.table('/home/koen/Data/SingleCell_Data/10X_Genomics/cluster_composition/653k_graphclust_clusters.csv_composition', quote='\'', comment.char='', col.names = c('N', 'ID'))
clus_comp_half<-clus_comp_half[!(clus_comp_half$ID=='Cluster'),]
clus_comp_half$ID <- as.numeric(clus_comp_half$ID)
clus_comp_half$percentage = clus_comp_half$N/sum(clus_comp_half$N)
clus_comp_half <- clus_comp_half[order(clus_comp_half$ID),]
clus_comp_half$percentage <- clus_comp_half$percentage*100


########## 1/8 of cells (n=163,277) ##########
clus_comp_eighth <- read.table('/home/koen/Data/SingleCell_Data/10X_Genomics/cluster_composition/163k_graphclust_clusters.csv_composition', quote='\'', comment.char='', col.names = c('N', 'ID'))
clus_comp_eighth<-clus_comp_eighth[!(clus_comp_eighth$ID=='Cluster'),]
clus_comp_eighth$ID <- as.numeric(clus_comp_eighth$ID)
clus_comp_eighth$percentage = clus_comp_eighth$N/sum(clus_comp_eighth$N)
clus_comp_eighth <- clus_comp_eighth[order(clus_comp_eighth$ID),]
clus_comp_eighth$percentage <- clus_comp_eighth$percentage*100


########## 1/12th of cells (n=108,851) ##########
clus_comp_twelfth <- read.table('/home/koen/Data/SingleCell_Data/10X_Genomics/cluster_composition/108k_graphclust_clusters.csv_composition', quote='\'', comment.char='', col.names = c('N', 'ID'))
clus_comp_twelfth<-clus_comp_twelfth[!(clus_comp_twelfth$ID=='Cluster'),]
clus_comp_twelfth$ID <- as.numeric(clus_comp_twelfth$ID)
clus_comp_twelfth$percentage = clus_comp_twelfth$N/sum(clus_comp_twelfth$N)
clus_comp_twelfth <- clus_comp_twelfth[order(clus_comp_twelfth$ID),]
clus_comp_twelfth$percentage <- clus_comp_twelfth$percentage*100


########## Plot cluster compositions ##########
pdf('cluster_compositions.pdf')
boxplot(percentage~ID, clus_comp_full, xlab = 'Cluster ID', ylab = 'Percentage of total', main = 'Cluster composition of full dataset (n=1,306,217)')
boxplot(percentage~ID, clus_comp_half, xlab = 'Cluster ID', ylab = 'Percentage of total', main = 'Cluster composition of 1/2 dataset (n=653,108)')
boxplot(percentage~ID, clus_comp_eighth, xlab = 'Cluster ID', ylab = 'Percentage of total', main = 'Cluster composition of 1/8 dataset (n=108,851)')
boxplot(percentage~ID, clus_comp_twelfth, xlab = 'Cluster ID', ylab = 'Percentage of total', main = 'Cluster composition of 1/12 dataset (n=108,851)')
boxplot(clus_comp_full$percentage, clus_comp_half$percentage, clus_comp_eighth$percentage, clus_comp_twelfth$percentage, xlab = 'Dataset', ylab = 'Percentage of total', main = 'Cluster composition of 10x Genomics subsampled datasets', names=c('Full', '1/2', '1/8', '1/12'), col=c('red', 'purple', 'cyan', 'pink'))
dev.off()
