#!/bin/bash
    #SBATCH -n 1
    #SBATCH -t hh:mm:ss
#SBATCH -p normal
    #SBATCH --mem=??G
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=K.J.Rademaker-2@umcutrecht.nl


    # ABOUT: Script to
    # REQUIRED:
# AUTHOR: Koen Rademaker
    # DATE:


# (0) File and variable organization
work_dir="$TMPDIR"/10x_genomics_qc_pp # Create work directory
mkdir ${work_dir}
out_dir="$TMPDIR"/10x_genomics_qc_pp/restructure_output # Create output directory
mkdir ${out_dir}

cp $HOME/Koen/SC_data/10x_Genomics/1M_neurons_filtered_gene_bc_matrices_h5.h5 ${work_dir} # Copy gene / cell matrix to work directory
cp $HOME/Koen/SC_data/10x_Genomics/analysis/clustering/graphclust/clusters.csv ${work_dir} # Copy clustering output to work directory
cp $HOME/Koen/SC_data/10x_Genomics/ensembl_v83_ensembl_mgisymbol_Mm.tab.gz ${work_dir}# Copy mouse mapping to work directory
cp $HOME/Koen/SC_data/10x_Genomics/ensembl_v82_Mm_Hs.tab.gz ${work_dir} # Copy mouse-to-human mapping to work directory
cd ${work_dir} # Move to work directory

${h5_file}=1M_neurons_filtered_gene_bc_matrices_h5.h5
${cluster_file}=clusters.csv
${mouse_mapping_file}=ensembl_v83_ensembl_mgisymbol_Mm.tab.gz
${mouse_human_mapping_file}=ensembl_v82_Mm_Hs.tab.gz

# (1) Execute Python program with the following arguments:
# - Gene / cell matrix - ${h5_file}
# - Clustering output - ${cluster_file}
# - Mouse mapping - ${mouse_mapping_file}
# - Mouse-to-human mapping - ${mouse_human_mapping_file}
    $HOME/Koen/Python-3.6.0/python ${script_path} ${h5_file} ${cluster_file} ${mouse_mapping_file} ${mouse_human_mapping_file}


# (3) File export and cleaning
    # DO
    # PLACE FILES INTO OUTPUT FOLDER WHERE NECESSARY
    cp ${out_dir} $HOME/Koen/SC_data/10x_Genomics/{NAME}
