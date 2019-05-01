#!/bin/bash

# GOAL:     Convert Karolinska Institue (KI) cluster-level mean gene expression data to DEPICT input
# INPUT:    KI data as .rda files (see README)
# AUTHOR:   Koen Rademaker
# DATE:     1 May 2019

########## Run R script to generate cluster-level mean expression files ##########
Rscript ki_to_cluster_mean_exp_matrix.R

########## Run Python script to restructure data for DEPICT analysis ##########
    # TODO
