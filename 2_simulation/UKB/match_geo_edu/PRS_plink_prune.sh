#!/bin/bash
#SBATCH --job-name=test
#SBATCH --mem=25G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:20:00

####################
/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--threads 5 \
--rm-dup exclude-all \
--mach-r2-filter \
--geno 0.001 \
--bfile /dcs04/nilanjan/data/zwang/family/simulation/pgs/data/mergedplink \
--score /dcs04/nilanjan/data/zwang/family/simulation/pgs/ea_score_pruned.txt 1 2 3 list-variants cols=+scoresums,-scoreavgs \
--out /dcs04/nilanjan/data/zwang/family/simulation/pgs/ea_score_parents_pruned

