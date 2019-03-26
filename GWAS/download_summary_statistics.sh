#!/bin/bash
# ABOUT: Script to automatically download, extract and organize summary statistics files.
# AUTHOR: Koen Rademaker
# DATE: 22 March 2019

# (0) Set up variables
scz_dir=scz_gwas
control_dir=control_gwas

# (1) Download and rename the files.
wget https://www.med.unc.edu/pgc/results-and-downloads/data-use-agreement-forms/copy2_of_mdd2018_ex23andme%20_data_download_agreement/resolveuid/a9451724b03741a4bf708c0014d7ff12 --no-check-certificate -O mdd_2018_ex23andme_sum_stats.txt.gz
wget https://ctg.cncr.nl/documents/p1651/SavageJansen_IntMeta_sumstats.zip --no-check-certificate -O intelligence_2018.zip
wget https://www.dropbox.com/s/ho58e9jmytmpaf8/GWAS_EA_excl23andMe.txt --no-check-certificate -O ea_2018_sum_stats.txt
wget http://cnsgenomics.com/data/yengo_et_al_2018_hmg/Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz --no-check-certificate -O he_2018_sum_stats.txt.gz
wget http://cnsgenomics.com/data/yengo_et_al_2018_hmg/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz --no-check-certificate -O bmi_2018_sum_stats.txt.gz
wget https://walters.psycm.cf.ac.uk/clozuk_pgc2.meta.sumstats.txt.gz --no-check-certificate -O scz_2018_sum_stats.txt.gz
wget https://walters.psycm.cf.ac.uk/clozuk_pgc2.meta.sumstats.info9.snplist.txt.gz --no-check-certificate -O scz_2018_hq_imputed_snps.txt.gz

# (2) Extract files
gunzip -v mdd_2018_ex23andme_sum_stats.txt.gz
unzip -p intelligence_2018.zip sumstats/SavageJansen_2018_intelligence_metaanalysis.txt > iq_2018_sum_stats.txt
gunzip -v he_2018_sum_stats.txt.gz
gunzip -v bmi_2018_sum_stats.txt.gz
gunzip -v scz_2018_sum_stats.txt.gz
gunzip -v scz_2018_hq_imputed_snps.txt.gz

# (3) Organize files
mkdir ${scz_dir}
mv scz_2018* ${scz_dir}
mkdir ${control_dir}
mv *.txt ${control_dir}
rm -rf intelligence_2018.zip
