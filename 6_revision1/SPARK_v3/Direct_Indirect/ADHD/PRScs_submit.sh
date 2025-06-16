#!/bin/bash
#SBATCH --job-name=family_child
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --time=23:30:00

python /dcs04/nilanjan/data/zwang/tools/PRScs/PRScs.py --ref_dir=/dcs04/nilanjan/data/zwang/phewas/LD/ldblk_ukbb_eur --bim_prefix=/dcs04/nilanjan/data/zwang/family/autism/bim_all_match1KG --sst_file=/dcs04/nilanjan/data/zwang/PRS_adhd/prscs_results/sum_stats.txt --n_gwas=128214 --out_dir=/dcs04/nilanjan/data/zwang/family/pleiotropy/redo/adhd/prscs_results/eur