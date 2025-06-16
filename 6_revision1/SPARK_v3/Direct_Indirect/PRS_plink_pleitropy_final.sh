#!/bin/bash
#SBATCH --job-name=test
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --array=1-10
#SBATCH --cpus-per-task=5
#SBATCH --time=23:30:00

####################
readarray -t a < /dcs04/nilanjan/data/zwang/family/pleiotropy/redo/pgs_trait.txt

/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--threads 5 \
--rm-dup exclude-all \
--mach-r2-filter \
--bfile /dcs05/ladd/NDEpi/data/Projects/InProgress/ziqiao/data/mergefinal  \
--read-freq /dcs05/ladd/NDEpi/data/Projects/InProgress/ziqiao/data/mergefinal.frq.counts \
--score /dcs04/nilanjan/data/zwang/family/pleiotropy/redo/${a[$(($SLURM_ARRAY_TASK_ID-1))]}_score.txt 1 2 3 list-variants cols=+scoresums,-scoreavgs \
--out /dcs04/nilanjan/data/zwang/family/pleiotropy/redo/${a[$(($SLURM_ARRAY_TASK_ID-1))]}/${a[$(($SLURM_ARRAY_TASK_ID-1))]}_score

