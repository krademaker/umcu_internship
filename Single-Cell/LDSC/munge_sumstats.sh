# GOAL:		Munge GWAS summary statistics for LDSC partitioned heritability
# REQUIRED:	GWAS summary statistics files, functioning LDSC environment (see README)
# AUTHOR:	Koen Rademaker
# DATE:		25 June 2019


########## Set up environment and variables ##########
ldsc_dir=./ldsc
sum_stats_dir=./Summary-Statistics
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2 -P Files/
bzip2 -d Files/w_hm3.snplist.bz2
cp Files/w_hm3.snplist ${ldsc_dir}
conda activate ldsc


########## Munge summary statistics ##########
# fix issue for missing items
python munge_sumstats.py --sumstats ~/Data/GWAS/original_sum_stats/clozuk_pgc2.meta.sumstats.txt.gz --merge-alleles ~/Downloads/tmp/w_hm3.snplist --N-cas 40675 --N-con 64643 --out SCZ

python ${ldsc_dir}/munge_sumstats.py --sumstats ${sum_stats_dir}/XXX --merge-alleles w_hm3.snplist --out XXX

python ${ldsc_dir}/munge_sumstats.py --sumstats ${sum_stats_dir}/XXX --merge-alleles w_hm3.snplist --out XXX

python ${ldsc_dir}/munge_sumstats.py --sumstats ${sum_stats_dir}/XXX --merge-alleles w_hm3.snplist --out XXX

python ${ldsc_dir}/munge_sumstats.py --sumstats ${sum_stats_dir}/XXX --merge-alleles w_hm3.snplist --out XXX

python ${ldsc_dir}/munge_sumstats.py --sumstats ${sum_stats_dir}/XXX --merge-alleles w_hm3.snplist --out XXX
# Munge for SCZ CLOZUK 2018
# Munge for BMI 2015
# Munge for MDD 2018
# Munge for EA 2018
# Munge for IQ 2018
# Munge for BC 2017
