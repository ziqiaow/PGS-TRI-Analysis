#!/bin/bash
#SBATCH --job-name=crc
#SBATCH --mem=50G
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1



/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440/plink/data_pfile/chr${SLURM_ARRAY_TASK_ID} \
--freq --out /dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440/plink/data_pfile/af_chr${SLURM_ARRAY_TASK_ID}

#https://www.biostars.org/p/306623/
#https://www.biostars.org/p/9495654/#9497482
#https://www.cog-genomics.org/plink/1.9/data