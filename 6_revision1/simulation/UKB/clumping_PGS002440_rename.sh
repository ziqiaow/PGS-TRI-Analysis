#!/bin/bash
#SBATCH --job-name=test
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --array=1-22
#SBATCH --cpus-per-task=6


####################
/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/data/chr${SLURM_ARRAY_TASK_ID} \
--keep /dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/ukb_eur_unrelated_ind.txt \
--clump /dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/sum_stats_chr${SLURM_ARRAY_TASK_ID}.txt \
--clump-p1 1 \
--clump-p2 1 \
--clump-r2 0.05 \
--clump-kb 250   \
--out /dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/snps_clumped_chr${SLURM_ARRAY_TASK_ID}

#--score-col-nums 3 \
#https://zzz.bwh.harvard.edu/plink/profile.shtml


#--extract /users/zwang4/GXE/ukbiobank/data/bc_313snp_test1.txt \