# TITLE:    MAGMA_Celltyping_10x_Genomics.R
# ABOUT:    Script to run MAGMA Celltyping analysis on 10x Genomics mouse brain scRNA-seq data
# INPUT:    mean_expr_path: Mean gene expression data for 10x Genomics cell types
#           scz_sum_stats: SCZ 2018 summary statistics, column order: SNP, CHR, BP, Freq.A1, A1, A2, OR, SE, P
#           mdd_sum_stats: MDD 2018 summary statistics, column order: SNP, CHR, BP, A1, A2, OR, SE, P
#           ea_sum_stats: Educational attainment 2018 summary statistics, column order: SNP, CHR, BP, A1, A2, EAF, Beta, SE, P
#           iq_sum_stats: Intelligence 2018 summary statistics, column order: SNP, CHR, POS, A1, A2, EAF_HRC, P, N
#           bmi_sum_stats: Body Mass Index 2015 summary statistics, column order: SNP, CHR, BP, A1, A2, MAF, B, SE, P, N
#           bc_sum_stats: Breast cancer 2017 summary statistics, column order: SNP, CHR, BP, A1, A2, MAF, BETA, SE, P
# AUTHOR:   Koen Rademaker
# DATE:     17 June 2019


########## Define custom functions ##########
custom_plot_celltype_associations <- function(ctAssocs,ctd,useSignificanceLine=TRUE,savePDF=TRUE,fileTag='',plotDendro=TRUE,gwas_title='',plotLegend=TRUE,figsDir=NA){
    # CHECK: THAT RESULTS FOR ONLY ONE GWAS WERE PROVIDED (for more than one use magma.tileplot.r)
    whichGWAS = unique(gsub('DOWN\\..*','DOWN',unique(ctAssocs[[1]]$results$GCOV_FILE)))
    if(length(whichGWAS)>1){stop('Only results for one GWAS at a tile should be provided to plot_celltype_association. For multiple GWAS, use magma.tileplot()')}
    magmaPaths = get.magma.paths(ctAssocs$gwas_sumstats_path,ctAssocs$upstream_kb,ctAssocs$downstream_kb)
    if(is.na(figsDir)){
        figsDir = magmaPaths$figs
    }
    # CHECK: WAS A TITLE PROVIDED FOR THE PLOT?
    if(gwas_title==''){gwas_title=whichGWAS}
    # CHECK: THAT A MINIMAL SET OF COLUMN HEADERS ARE INCLUDED IN THE RESULTS TABLE
    requiredHeaders = c('Celltype','P','log10p','Method','EnrichmentMode','CONTROL','CONTROL_label')
    # Is the analysis top10%, linear or merged?
    print(ctAssocs[[1]]$results)
    print(unique(ctAssocs[[1]]$results$EnrichmentMode))
    if(length(unique(ctAssocs[[1]]$results$EnrichmentMode))==1){
        if(unique(ctAssocs[[1]]$results$EnrichmentMode)=='Linear'){
            analysisType = 'Linear'
        }else{
            analysisType = 'TopDecile'
        }
    }else{  analysisType = 'Merged' }
    # Generate the plots (for each annotation level seperately)
    library(cowplot)
    figures = list()
    for(annotLevel in 1:sum(names(ctAssocs)=='')){
        # SET: NEW COLUMN COMBINING METHODS or ENRICHMENTTYPES
        ctAssocs[[annotLevel]]$results$FullMethod = sprintf('%s %s',ctAssocs[[annotLevel]]$results$Method,ctAssocs[[annotLevel]]$results$EnrichmentMode)
        if(plotDendro==TRUE){
            # Order cells by dendrogram
            ctdDendro = get.ctd.dendro(ctd,annotLevel=annotLevel)
            ctAssocs[[annotLevel]]$results$Celltype <- factor(ctAssocs[[annotLevel]]$results$Celltype, levels=gsub(' |\\(|\\)','\\.',ctdDendro$ordered_cells))
        }
        a2 <- ggplot(ctAssocs[[annotLevel]]$results, aes(x = factor(Celltype), y = log10p, fill=FullMethod)) + scale_y_reverse()+geom_bar(stat = 'identity',position='dodge') + coord_flip() + ylab(expression('log'[10]*'(pvalue)')) + xlab('')
        a2 <- a2 + theme(legend.position = c(0.5, 0.8)) + ggtitle(gwas_title) + theme(legend.title=element_blank())
        if(plotLegend==FALSE){    a2 = a2 + theme(legend.position='none') }
        if(useSignificanceLine){  a2 = a2+geom_hline(yintercept=log(as.numeric(0.05/ctAssocs$total_baseline_tests_performed),10),colour='black')    }
        theFig = a2
        # If the results come from a BASELINE analysis... 
        if(length(unique(ctAssocs[[1]]$results$CONTROL))==1){
            if(plotDendro==TRUE){     theFig = grid.arrange(a2,ctdDendro$dendroPlot,ncol=2,widths=c(0.8,0.2))        }
            if(savePDF){
                fName = sprintf('%s/%s.%sUP.%sDOWN.annotLevel%s.Baseline.%s.%s.pdf',figsDir,magmaPaths$gwasFileName,ctAssocs$upstream_kb,ctAssocs$downstream_kb,annotLevel,fileTag,analysisType)
                print('here')
                pdf(file=fName,width=10,height=1+2*(dim(ctAssocs[[annotLevel]]$results)[1]/10))
                print(grid.arrange(a2,ctdDendro$dendroPlot,ncol=2,widths=c(0.8,0.2)))
                dev.off()
            }else{print(theFig)}            
            # IF THE RESULTS COME FROM A CONDITIONAL ANALYSIS    
        }else{
            theFig = theFig + facet_wrap(~CONTROL_label)
            if(savePDF){
                fName = sprintf('%s/%s.%sUP.%sDOWN.annotLevel%s.ConditionalFacets.%s.%s.pdf',figsDir,magmaPaths$gwasFileName,ctAssocs$upstream_kb,ctAssocs$downstream_kb,annotLevel,fileTag,analysisType)
                pdf(file=fName,width=25,height=1+2*(dim(ctAssocs[[annotLevel]]$results)[1]/30))
                print(theFig)
                dev.off()
            }else{print(theFig)}            
        }
        figures[[length(figures)+1]] = theFig
    }
    return(figures)
}


########## Load required libraries ##########
library(MAGMA.Celltyping)
library(One2One)


########## Set variables ##########
genome_ref_dir = '~/Git/umcu_internship/Single-Cell/MAGMA-Celltyping/g1000_eur'
if(!file.exists(sprintf('%s/g1000_eur.bed', genome_ref_dir))){
    download.file('https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip',destfile=sprintf('%s.zip', genome_ref_dir))
    unzip(sprintf('%s.zip',genome_ref_dir), exdir=genome_ref_dir)
}
genome_ref_path = sprintf('%s/g1000_eur', genome_ref_dir)
sum_stats_dir = '~/Git/umcu_internship/Single-Cell/MAGMA-Celltyping/Summary-Statistics/'
suppl_data_dir = '~/Data/MAGMA-Celltyping-GWAS/supplementary_data/'

mean_expr_path <- '~/Git/umcu_internship/Single-Cell/data/10x_Genomics_16_cell_types_mean_expression.tsv'
scz_sum_stats = paste0(sum_stats_dir, 'clozuk_pgc2.meta.sumstats_MAGMA.txt')
mdd_sum_stats = paste0(sum_stats_dir, 'MDD2018_ex23andMe_MAGMA.txt')
ea_sum_stats = paste0(sum_stats_dir, 'EA_excl23andMe_MAGMA.txt')
iq_sum_stats = paste0(sum_stats_dir, 'Intelligence_2018_MAGMA.txt')
bmi_sum_stats = paste0(sum_stats_dir, 'BMI_2015_MAGMA.txt')
bc_sum_stats = paste0(sum_stats_dir, 'BreastCancer_2017_MAGMA.txt')

formal_cell_type_names <- c('Glutamatergic neurons', 'Neuroblasts 1', 'Astrocytes 1', 'Neuroblasts 2', 'Intermediate progenitors', 'Enteric glial cells', 'Interneurons 1', 'Neurons', 'Interneurons 2', 'Vascular endothelial cells', 'Enteric neurons', 'Cajal-Retzius cells', 'Oligodendrocytes', 'Astrocytes 2', 'Endothelial cells', 'Microglia')


########## Load cell type-level mean expression data & calculate specificity metric S(g,c) ##########
mean_expr_data <- read.delim(mean_expr_path, row.names = 1)
specificity_data <- mean_expr_data/rowSums(mean_expr_data)


########## Restructure ctd object for MAGMA cell type association analysis ##########
ctd_10X <- MAGMA.Celltyping::ctd_allKI
ctd_10X[[1]]$mean_exp <- mean_expr_data
ctd_10X[[1]]$specificity <- specificity_data


########## Calculate specificity quantiles in 41 bins (least-to-most specific expression) ##########
ctd_10X = prepare.quantile.groups(ctd_10X, specificity_species = 'mouse', numberOfBins = 41)
ctd_10X[[2]] <- NULL

########## Optional plots of specificity ##########
ggplot(ctd_10X[[1]]$specificity, aes(x = Interneurons_1)) + geom_histogram(bins = 100)


#################### Analysis for Pardiñas et al. 2018 schizophrenia GWAS ####################
# (1) Map SNPs to genes
scz_mapped_genes = map.snps.to.genes(scz_sum_stats, genome_ref_path = genome_ref_path, N = 105318)
# (2) Perform linear cell type association analysis
scz_linear_assoc = calculate_celltype_associations(ctd_10X, scz_sum_stats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Linear')
# (3) Perform top 10% cell type association analysis
scz_top_10_assoc = calculate_celltype_associations(ctd_10X, scz_sum_stats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Top 10%')
# (4) Update cell type names
scz_linear_assoc[[1]]$results$Celltype = formal_cell_type_names
scz_top_10_assoc[[1]]$results$Celltype = formal_cell_type_names
# (5) Save analysis results as files
write.table(scz_linear_assoc[[1]]$results, file = paste0(suppl_data_dir, '10x_lvl1_linear_scz_2018.tsv'), sep = '\t', row.names = FALSE)
write.table(scz_top_10_assoc[[1]]$results, file = paste0(suppl_data_dir, '10x_lvl1_top10_scz_2018.tsv'), sep = '\t', row.names = FALSE)
# (6) Plot analysis results
scz_linear_figures = custom_plot_celltype_associations(scz_linear_assoc, ctd = ctd_10X, useSignificanceLine = FALSE, plotLegend = FALSE, plotDendro = FALSE, gwas_title = 'Linear cell type association analysis - SCZ 2018', savePDF = FALSE)
scz_top_10_figures = custom_plot_celltype_associations(scz_top_10_assoc, ctd = ctd_10X, useSignificanceLine = FALSE, plotLegend = FALSE, plotDendro = FALSE, gwas_title = 'Top 10% cell type association analysis - SCZ 2018', savePDF = FALSE)


#################### Analysis for Wray et al. 2018 Major Depressive Disorder (MDD) GWAS ####################
# (1) Map SNPs to genes
mdd_mapped_genes = map.snps.to.genes(mdd_sum_stats, genome_ref_path = genome_ref_path, N = 480359)
# (2) Perform linear cell type association analysis
mdd_linear_assoc = calculate_celltype_associations(ctd_10X, mdd_sum_stats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Linear')
# (3) Perform top 10% cell type association analysis
mdd_top_10_assoc = calculate_celltype_associations(ctd_10X, mdd_sum_stats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Top 10%')
# (4) Update cell type names
mdd_linear_assoc[[1]]$results$Celltype = formal_cell_type_names
mdd_top_10_assoc[[1]]$results$Celltype = formal_cell_type_names
# (5) Save analysis results as files
write.table(mdd_linear_assoc[[1]]$results, file = paste0(suppl_data_dir, '10x_lvl1_linear_mdd_2018.tsv'), sep = '\t', row.names = FALSE)
write.table(mdd_top_10_assoc[[1]]$results, file = paste0(suppl_data_dir, '10x_lvl1_top10_mdd_2018.tsv'), sep = '\t', row.names = FALSE)
# (6) Plot analysis results
mdd_linear_figures = custom_plot_celltype_associations(mdd_linear_assoc, ctd = ctd_10X, useSignificanceLine = FALSE, plotLegend = FALSE, plotDendro = FALSE, gwas_title = 'Linear cell type association analysis - MDD 2018', savePDF = FALSE)
mdd_top_10_figures = custom_plot_celltype_associations(mdd_top_10_assoc, ctd = ctd_10X, useSignificanceLine = FALSE, plotLegend = FALSE, plotDendro = FALSE, gwas_title = 'Top 10% cell type association analysis - MDD 2018', savePDF = FALSE)


#################### Analysis for Lee et al. 2018 Educational Attainment (EA) GWAS ####################
# (1) Map SNPs to genes
ea_mapped_genes = map.snps.to.genes(ea_sum_stats, genome_ref_path = genome_ref_path, N = 1318881)
# (2) Perform linear cell type association analysis
ea_linear_assoc = calculate_celltype_associations(ctd_10X, ea_sum_stats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Linear')
# (3) Perform top 10% cell type association analysis
ea_top_10_assoc = calculate_celltype_associations(ctd_10X, ea_sum_stats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Top 10%')
# (4) Update cell type names
ea_linear_assoc[[1]]$results$Celltype = formal_cell_type_names
ea_top_10_assoc[[1]]$results$Celltype = formal_cell_type_names
# (5) Save analysis results as files
write.table(ea_linear_assoc[[1]]$results, file = paste0(suppl_data_dir, '10x_lvl1_linear_ea_2018.tsv'), sep = '\t', row.names = FALSE)
write.table(ea_top_10_assoc[[1]]$results, file = paste0(suppl_data_dir, '10x_lvl1_top10_ea_2018.tsv'), sep = '\t', row.names = FALSE)
# (6) Plot analysis results
ea_linear_figures = custom_plot_celltype_associations(ea_linear_assoc, ctd = ctd_10X, useSignificanceLine = FALSE, plotLegend = FALSE, plotDendro = FALSE, gwas_title = 'Linear cell type association analysis - EA 2018', savePDF = FALSE)
ea_top_10_figures = custom_plot_celltype_associations(ea_top_10_assoc, ctd = ctd_10X, useSignificanceLine = FALSE, plotLegend = FALSE, plotDendro = FALSE, gwas_title = 'Top 10% cell type association analysis - EA 2018', savePDF = FALSE)


#################### Analysis for Savage et al. 2018 intelligence GWAS ####################
# (1) Map SNPs to genes
iq_mapped_genes = map.snps.to.genes(iq_sum_stats, genome_ref_path = genome_ref_path, N = 269897)
# (2) Perform linear cell type association analysis
iq_linear_assoc = calculate_celltype_associations(ctd_10X, iq_sum_stats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Linear')
# (3) Perform top 10% cell type association analysis
iq_top_10_assoc = calculate_celltype_associations(ctd_10X, iq_sum_stats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Top 10%')
# (4) Update cell type names
iq_linear_assoc[[1]]$results$Celltype = formal_cell_type_names
iq_top_10_assoc[[1]]$results$Celltype = formal_cell_type_names
# (5) Save analysis results as files
write.table(iq_linear_assoc[[1]]$results, file = paste0(suppl_data_dir, '10x_lvl1_linear_iq_2018.tsv'), sep = '\t', row.names = FALSE)
write.table(iq_top_10_assoc[[1]]$results, file = paste0(suppl_data_dir, '10x_lvl1_top10_iq_2018.tsv'), sep = '\t', row.names = FALSE)
# (6) Plot analysis results
iq_linear_figures = custom_plot_celltype_associations(iq_linear_assoc, ctd = ctd_10X, useSignificanceLine = FALSE, plotLegend = FALSE, plotDendro = FALSE, gwas_title = 'Linear cell type association analysis - Intelligence 2018', savePDF = FALSE)
iq_top_10_figures = custom_plot_celltype_associations(iq_top_10_assoc, ctd = ctd_10X, useSignificanceLine = FALSE, plotLegend = FALSE, plotDendro = FALSE, gwas_title = 'Top 10% cell type association analysis - Intelligence 2018', savePDF = FALSE)


#################### Analysis for Locke et al. 2015 Body Mass Index (BMI) GWAS ####################
# (1) Map SNPs to genes
bmi_mapped_genes = map.snps.to.genes(bmi_sum_stats, genome_ref_path = genome_ref_path)
# (2) Perform linear cell type association analysis
bmi_linear_assoc = calculate_celltype_associations(ctd_10X, bmi_sum_stats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Linear')
# (3) Perform top 10% cell type association analysis
bmi_top_10_assoc = calculate_celltype_associations(ctd_10X, bmi_sum_stats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Top 10%')
# (4) Update cell type names
bmi_linear_assoc[[1]]$results$Celltype = formal_cell_type_names
bmi_top_10_assoc[[1]]$results$Celltype = formal_cell_type_names
# (5) Save analysis results as files
write.table(bmi_linear_assoc[[1]]$results, file = paste0(suppl_data_dir, '10x_lvl1_linear_bmi_2015.tsv'), sep = '\t', row.names = FALSE)
write.table(bmi_top_10_assoc[[1]]$results, file = paste0(suppl_data_dir, '10x_lvl1_top10_bmi_2015.tsv'), sep = '\t', row.names = FALSE)
# (6) Plot analysis results
bmi_linear_figures = custom_plot_celltype_associations(bmi_linear_assoc, ctd = ctd_10X, useSignificanceLine = FALSE, plotLegend = FALSE, plotDendro = FALSE, gwas_title = 'Linear cell type association analysis - BMI 2015', savePDF = FALSE)
bmi_top_10_figures = custom_plot_celltype_associations(bmi_top_10_assoc, ctd = ctd_10X, useSignificanceLine = FALSE, plotLegend = FALSE, plotDendro = FALSE, gwas_title = 'Top 10% cell type association analysis - BMI 2015', savePDF = FALSE)


#################### Analysis for Michailidou et al. 2017 breast cancer GWAS ####################
# (1) Map SNPs to genes
bc_mapped_genes = map.snps.to.genes(bc_sum_stats, genome_ref_path = genome_ref_path, N = 256123)
# (2) Perform linear cell type association analysis
bc_linear_assoc = calculate_celltype_associations(ctd_10X, bc_sum_stats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Linear')
# (3) Perform top 10% cell type association analysis
bc_top_10_assoc = calculate_celltype_associations(ctd_10X, bc_sum_stats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Top 10%')
# (4) Update cell type names
bc_linear_assoc[[1]]$results$Celltype = formal_cell_type_names
bc_top_10_assoc[[1]]$results$Celltype = formal_cell_type_names
# (5) Save analysis results as files
write.table(bc_linear_assoc[[1]]$results, file = paste0(suppl_data_dir, '10x_lvl1_linear_bc_2017.tsv'), sep = '\t', row.names = FALSE)
write.table(bc_top_10_assoc[[1]]$results, file = paste0(suppl_data_dir, '10x_lvl1_top10_bc_2017.tsv'), sep = '\t', row.names = FALSE)
# (6) Plot analysis results
bc_linear_figures = custom_plot_celltype_associations(bc_linear_assoc, ctd = ctd_10X, useSignificanceLine = FALSE, plotLegend = FALSE, plotDendro = FALSE, gwas_title = 'Linear cell type association analysis - Breast cancer 2017', savePDF = FALSE)
bc_top_10_figures = custom_plot_celltype_associations(bc_top_10_assoc, ctd = ctd_10X, useSignificanceLine = FALSE, plotLegend = FALSE, plotDendro = FALSE, gwas_title = 'Top 10% cell type association analysis - Breast cancer 2017', savePDF = FALSE)