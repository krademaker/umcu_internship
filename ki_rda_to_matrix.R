# TITLE:    ki_to_cluster_mean_exp_matrix.R
# ABOUT:    Script to save cluster-level (n=49 & n=149) mean expression data to files for further analysis.
# INPUT:    ctd_allKI.rda: RData object containing specificity and mean expression data
# OUTPUT:   lvl_1_mean_exp_matrix.txt: Matrix with mean expression data for level 1 clustering (n=49)
# OUTPUT:   lvl_2_mean_exp_matrix.txt: Matrix with mean expression data for level 2 clustering (n=149)
# AUTHOR:   Koen Rademaker
# DATE:     30 April 2019

########## Load RData object ##########
load('/home/koen/Git/umcu_internship/Single-Cell/ctd_allKI.rda')

########## Select level 1 data and save to file ##########
level_1_mean_expression_matrix <- ctd[[1]]$mean_exp
write.table(level_1_mean_expression_matrix,
            file = '/home/koen/Data/lvl_1_mean_exp_matrix.txt',
            row.names=T,
            col.names=NA,
            sep = '\t')

########## Select level 2 data and save to file ##########
level_2_mean_expression_matrix <- ctd[[2]]$mean_exp
write.table(level_2_mean_expression_matrix,
            file = '/home/koen/Data/lvl_2_mean_exp_matrix.txt',
            row.names=T,
            col.names=NA,
            sep = '\t')