#!/bin/bash

# GOAL:     Convert Karolinska Institue (KI) cluster-level mean gene expression data to DEPICT input
# INPUT:    KI data as .rda files (see README)
# AUTHOR:   Koen Rademaker
# DATE:     5 May 2019

script_path=~/Git/umcu_internship/Single-Cell/KI-Data-DEPICT
data_path=~/Data/SingleCell_Data/Skene-2018

########## Run R script to generate cluster-level mean expression files ##########
Rscript ${script_path}/ki_rda_to_matrix.R

########## Run Python script to restructure data for DEPICT analysis ##########
python ${script_path}/ki_matrix_to_depict.py ${data_path}/KI_mean_lvl1_exp.tsv ${data_path}/KI_lvl1_standardized.tsv
python ${script_path}/ki_matrix_to_depict.py ${data_path}/KI_mean_lvl2_exp.tsv ${data_path}/KI_lvl2_standardized.tsv
