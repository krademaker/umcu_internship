# GOAL:		Munge GWAS summary statistics for LDSC partitioned heritability
# REQUIRED:	GWAS summary statistics files, functioning LDSC environment (see README)
# AUTHOR:	Koen Rademaker
# DATE:		27 June 2019


########## Set up environment and variables ##########
ldsc_dir=./ldsc
sum_stats_dir=./Summary-Statistics
data_dir=./Files
conda activate ldsc


########## Munge summary statistics ##########

# Download list of HapMap3 SNPs
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2 -P ${data_dir}
bzip2 -d ${data_dir}/w_hm3.snplist.bz2


# Download & pre-process schizophrenia summary statistics
wget https://walters.psycm.cf.ac.uk/clozuk_pgc2.meta.sumstats.txt.gz --no-check-certificate -P ${sum_stats_dir}
gunzip ${sum_stats_dir}/clozuk_pgc2.meta.sumstats.txt.gz
head -n 1 ${sum_stats_dir}/clozuk_pgc2.meta.sumstats.txt > ${data_dir}/tmp_clozuk_header
grep 'rs' ${sum_stats_dir}/clozuk_pgc2.meta.sumstats.txt > ${sum_stats_dir}/clozuk_pgc2.meta.sumstats.only_rs.txt
awk -v n=1 '{s = gensub("^(([^:]*:){"n"}).*$", "\\1", 1); sub(".$","",s); print $0, s}' ${sum_stats_dir}/clozuk_pgc2.meta.sumstats.only_rs.txt > ${sum_stats_dir}/clozuk_pgc2.meta.sumstats.only_rs.rs_col.txt
sed -i 's/\r//' ${sum_stats_dir}/clozuk_pgc2.meta.sumstats.only_rs.rs_col.txt
sed -i 's/\r//' ${data_dir}/tmp_clozuk_header
sed 's/$/\tSNP_RS/' ${data_dir}/tmp_clozuk_header > ${data_dir}/tmp_clozuk_header_2
cat ${data_dir}/tmp_clozuk_header_2 ${sum_stats_dir}/clozuk_pgc2.meta.sumstats.only_rs.rs_col.txt > ${sum_stats_dir}/clozuk_pgc2.meta.sumstats.only_rs.rs_col.header.txt
# Munge schizophrenia summary statistics
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sum_stats_dir}/clozuk_pgc2.meta.sumstats.only_rs.rs_col.header.txt --out SCZ --N 105318 --merge-alleles ${data_dir}/w_hm3.snplist --a1-inc --snp SNP_RS --a1 A1 --a2 A2 --p P --frq Freq.A1


# Pre-process MDD summary statistics
gunzip ${sum_stats_dir}/MDD2018_ex23andMe.gz
# Munge MDD summary statistics
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sum_stats_dir}/MDD2018_ex23andMe --out MDD --N 173005 --merge-alleles ${data_dir}/w_hm3.snplist --snp SNP --a1 A1 --a2 A2 --p P --info INFO --frq FRQ_U_113154


# Pre-process intelligence summary statistics
unzip ${sum_stats_dir}/SavageJansen_IntMeta_sumstats.zip
# Munge intelligence summary statistics
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sum_stats_dir}/sumstats/SavageJansen_2018_intelligence_metaanalysis.txt --out Int --n-col N_analyzed --merge-alleles ${data_dir}/w_hm3.snplist --snp SNP --a1 A1 --a2 A2 --frq EAF_HRC --info minINFO


# Munge educational attainment summary statistics
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sum_stats_dir}/GWAS_EA_excl23andMe.txt --out EA --N 1318881 --merge-alleles ${data_dir}/w_hm3.snplist --snp MarkerName --a1 A1 --a2 A2 --frq EAF


# Pre-process breast cancer summary statistics
gunzip ${sum_stats_dir}/oncoarray_bcac_public_release_oct17\ \(1\).txt.gz
mv ${sum_stats_dir}/oncoarray_bcac_public_release_oct17\ \(1\).txt ${sum_stats_dir}/oncoarray_bcac_public_release_oct17.txt
cut -f2-10 ${sum_stats_dir}/oncoarray_bcac_public_release_oct17.txt > ${sum_stats_dir}/oncoarray_bcac_public_release_oct17.main.txt
head -n 1 ${sum_stats_dir}/oncoarray_bcac_public_release_oct17.main.txt > ${data_dir}/tmp_oncoarray_bcac_header
grep 'rs' ${sum_stats_dir}/oncoarray_bcac_public_release_oct17.main.txt > ${sum_stats_dir}/oncoarray_bcac_public_release_oct17.main.only_rs.txt
awk -v n=1 '{s = gensub("^(([^:]*:){"n"}).*$", "\\1", 1); sub(".$","",s); print $0, s}' ${sum_stats_dir}/oncoarray_bcac_public_release_oct17.main.only_rs.txt > ${sum_stats_dir}/oncoarray_bcac_public_release_oct17.main.only_rs.rs_col.txt
sed -i 's/\r//' ${sum_stats_dir}/oncoarray_bcac_public_release_oct17.main.only_rs.rs_col.txt
sed -i 's/\r//' ${data_dir}/tmp_oncoarray_bcac_header
sed 's/$/\tSNP_RS/' ${data_dir}/tmp_oncoarray_bcac_header > ${data_dir}/tmp_oncoarray_bcac_header_2
cat ${data_dir}/tmp_oncoarray_bcac_header_2 ${sum_stats_dir}/oncoarray_bcac_public_release_oct17.main.only_rs.rs_col.txt > ${sum_stats_dir}/oncoarray_bcac_public_release_oct17.main.only_rs.rs_col.header.txt
# Munge breast cancer summary statistics
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sum_stats_dir}/oncoarray_bcac_public_release_oct17.main.only_rs.rs_col.header.txt --out BreastCancer --N 256123 --merge-alleles ${data_dir}/w_hm3.snplist --a1-inc --snp SNP_RS --a1 a0 --a2 a1 --p bcac_onco_icogs_gwas_P1df --frq bcac_onco_icogs_gwas_eaf_controls


# Download & pre-process BMI summary statistics
wget https://portals.broadinstitute.org/collaboration/giant/images/f/f0/All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz -P ${sum_stats_dir}
gunzip ${sum_stats_dir}/All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz
# Munge BMI summary statistics
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sum_stats_dir}/All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq --out BMI --merge-alleles ${data_dir}/w_hm3.snplist --snp SNP --a1 A1 --a2 A2 --frq Freq1.Hapmap -p p
