# TITLE:    cluster_composition.R
# ABOUT:    Script to determine and visualize cluster composition for the full or subsampled 10x Genomics mouse brain cell dataset.
# INPUT:    [N]_graphcluster_clusters.csv_composition: File(s) detailing number of cells per cluster (N) and cluster ID (ID)
# AUTHOR:   Koen Rademaker
# DATE:     14 May 2019


########## All cells (n=1,306,127) ##########
clust_comp_1m <- read.table('/home/koen/Data/SingleCell_Data/10X_Genomics/cluster_composition/graphclust_clusters.csv_composition', quote='\'', comment.char='', col.names = c('N', 'ID'))
clust_comp_1m <- clust_comp_1m[!(clust_comp_1m$ID=='Cluster'),]
clust_comp_1m$ID <- as.numeric(clust_comp_1m$ID)
clust_comp_1m$percentage = clust_comp_1m$N/sum(clust_comp_1m$N)
clust_comp_1m <- clust_comp_1m[order(clust_comp_1m$ID),]
clust_comp_1m$percentage <- clust_comp_1m$percentage*100


########## 1/2 of cells (n=653,064) ##########
clust_comp_653k <- read.table('/home/koen/Data/SingleCell_Data/10X_Genomics/cluster_composition/653k_graphclust_clusters.csv_composition', quote='\'', comment.char='', col.names = c('N', 'ID'))
clust_comp_653k<-clust_comp_653k[!(clust_comp_653k$ID=='Cluster'),]
clust_comp_653k$ID <- as.numeric(clust_comp_653k$ID)
clust_comp_653k$percentage = clust_comp_653k$N/sum(clust_comp_653k$N)
clust_comp_653k <- clust_comp_653k[order(clust_comp_653k$ID),]
clust_comp_653k$percentage <- clust_comp_653k$percentage*100


########## 1/4 of cells (n=326,531) ##########
clust_comp_326k <- read.table('/home/koen/Data/SingleCell_Data/10X_Genomics/cluster_composition/326k_graphclust_clusters.csv_composition', quote='\'', comment.char='', col.names = c('N', 'ID'))
clust_comp_326k<-clust_comp_326k[!(clust_comp_326k$ID=='Cluster'),]
clust_comp_326k$ID <- as.numeric(clust_comp_326k$ID)
clust_comp_326k$percentage = clust_comp_326k$N/sum(clust_comp_326k$N)
clust_comp_326k <- clust_comp_326k[order(clust_comp_326k$ID),]
clust_comp_326k$percentage <- clust_comp_326k$percentage*100


########## 1/8 of cells (n=163,265) ##########
clust_comp_163k <- read.table('/home/koen/Data/SingleCell_Data/10X_Genomics/cluster_composition/163k_graphclust_clusters.csv_composition', quote='\'', comment.char='', col.names = c('N', 'ID'))
clust_comp_163k<-clust_comp_163k[!(clust_comp_163k$ID=='Cluster'),]
clust_comp_163k$ID <- as.numeric(clust_comp_163k$ID)
clust_comp_163k$percentage = clust_comp_163k$N/sum(clust_comp_163k$N)
clust_comp_163k <- clust_comp_163k[order(clust_comp_163k$ID),]
clust_comp_163k$percentage <- clust_comp_163k$percentage*100


########## 1/12th of cells (n=108,844) ##########
clust_comp_108k <- read.table('/home/koen/Data/SingleCell_Data/10X_Genomics/cluster_composition/108k_graphclust_clusters.csv_composition', quote='\'', comment.char='', col.names = c('N', 'ID'))
clust_comp_108k<-clust_comp_108k[!(clust_comp_108k$ID=='Cluster'),]
clust_comp_108k$ID <- as.numeric(clust_comp_108k$ID)
clust_comp_108k$percentage = clust_comp_108k$N/sum(clust_comp_108k$N)
clust_comp_108k <- clust_comp_108k[order(clust_comp_108k$ID),]
clust_comp_108k$percentage <- clust_comp_108k$percentage*100


########## Plot cluster compositions ##########
pdf('cluster_compositions.pdf')
boxplot(percentage~ID, clust_comp_1m, xlab = 'Cluster ID', ylab = 'Percentage of total', main = 'Cluster composition of full dataset (n=1,306,127)', ylim = c(0.0, 10.0), las = 2, cex.axis=0.5)
boxplot(percentage~ID, clust_comp_653k, xlab = 'Cluster ID', ylab = 'Percentage of total', main = 'Cluster composition of 1/2 dataset (n=653,064)', ylim = c(0.0, 10.0), las = 2, cex.axis=0.5)
boxplot(percentage~ID, clust_comp_326k, xlab = 'Cluster ID', ylab = 'Percentage of total', main = 'Cluster composition of 1/4 dataset (n=326,531)', ylim = c(0.0, 10.0), las = 2, cex.axis=0.5)
boxplot(percentage~ID, clust_comp_163k, xlab = 'Cluster ID', ylab = 'Percentage of total', main = 'Cluster composition of 1/8 dataset (n=163,265)', ylim = c(0.0, 10.0), las = 2, cex.axis=0.5)
boxplot(percentage~ID, clust_comp_108k, xlab = 'Cluster ID', ylab = 'Percentage of total', main = 'Cluster composition of 1/12 dataset (n=108,844)', ylim = c(0.0, 10.0), las = 2, cex.axis=0.5)
boxplot(clust_comp_1m$percentage, clust_comp_653k$percentage, clust_comp_326k$percentage, clust_comp_163k$percentage, clust_comp_108k$percentage, xlab = 'Dataset', ylab = 'Percentage of total', main = 'Cluster composition of 10x Genomics subsampled datasets', names=c('Full', '1/2', '1/4', '1/8', '1/12'), col=c('red', 'purple', 'yellow', 'cyan', 'pink'))
dev.off()
