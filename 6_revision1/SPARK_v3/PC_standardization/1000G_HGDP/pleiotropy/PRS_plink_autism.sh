#!/bin/bash
#SBATCH --job-name=asd
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --array=1-22
#SBATCH --cpus-per-task=5


####################
/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--threads 5 \
--mach-r2-filter \
--rm-dup exclude-all \
--bfile /dcs04/nilanjan/data/ydun/HGDP+1000G/GRCh38/original_files/chr${SLURM_ARRAY_TASK_ID} \
--remove /dcs04/nilanjan/data/ydun/HGDP+1000G/GRCh38/related_low_quality.id \
--score /dcs04/nilanjan/data/zwang/family/autism/plink/pgs/redo/autism_score.txt 1 2 3 list-variants cols=+scoresums,-scoreavgs \
--out /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/pgs/score_chr$SLURM_ARRAY_TASK_ID  
