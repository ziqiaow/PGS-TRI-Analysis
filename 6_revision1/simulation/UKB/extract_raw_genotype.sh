#!/bin/bash
#SBATCH --job-name=crc
#SBATCH --mem=50G
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/data/chr${SLURM_ARRAY_TASK_ID} \
--extract /dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/clumped_snps_chr${SLURM_ARRAY_TASK_ID}.txt \
--recode A --out /dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/data/chr${SLURM_ARRAY_TASK_ID}

