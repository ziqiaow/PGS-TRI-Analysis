#!/bin/bash
#SBATCH --job-name=GTExv8.ALL.Brain_Hippocampus
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --array=1-3840
#SBATCH --cpus-per-task=1
#SBATCH --time=0-23:00:00

####################
r=GTExv8.ALL.Brain_Hippocampus

readarray -t a < /dcs04/nilanjan/data/zwang/family/gtex/data/weights_clean/${r}.txt

/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--threads 1 \
--rm-dup exclude-all \
--mach-r2-filter \
--bfile /dcs05/ladd/NDEpi/data/Projects/InProgress/ziqiao/data/mergefinal  \
--extract /dcs04/nilanjan/data/zwang/family/gtex/data/weights_clean/${r}/${a[$(($SLURM_ARRAY_TASK_ID-1))]}.txt \
--score /dcs04/nilanjan/data/zwang/family/gtex/data/weights_clean/${r}/${a[$(($SLURM_ARRAY_TASK_ID-1))]}.txt 1 2 3 list-variants cols=+scoresums,-scoreavgs \
--out /dcs04/nilanjan/data/zwang/family/gtex/data/scores/${r}/${a[$(($SLURM_ARRAY_TASK_ID-1))]}
