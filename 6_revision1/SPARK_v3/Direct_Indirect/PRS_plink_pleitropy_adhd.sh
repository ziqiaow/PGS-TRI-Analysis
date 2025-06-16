#!/bin/bash
#SBATCH --job-name=test
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=23:30:00

####################
/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--threads 1 \
--rm-dup exclude-all \
--mach-r2-filter \
--extract /dcs04/nilanjan/data/zwang/family/pleiotropy/redo/adhd_score.txt \
--bfile /dcs05/ladd/NDEpi/data/Projects/InProgress/ziqiao/data/mergefinal  \
--score /dcs04/nilanjan/data/zwang/family/pleiotropy/redo/adhd_score.txt 1 2 3 list-variants cols=+scoresums,-scoreavgs \
--out /dcs04/nilanjan/data/zwang/family/pleiotropy/redo/adhd/adhd_score


