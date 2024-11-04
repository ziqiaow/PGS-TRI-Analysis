#!/bin/bash
#SBATCH --job-name=crc
#SBATCH --mem=15G
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:30:00

/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--pfile /dcl01/arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr${SLURM_ARRAY_TASK_ID} \
--keep /dcs04/nilanjan/data/zwang/family/simulation/match_geo/plink/keep_parents.txt \
--extract range /dcs04/nilanjan/data/zwang/family/simulation/pgs/snpfile.txt \
--make-bed --out /dcs04/nilanjan/data/zwang/family/simulation/match_geo/plink/data/chr${SLURM_ARRAY_TASK_ID}

