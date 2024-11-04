#!/bin/bash
#SBATCH --job-name=prune
#SBATCH --mem=15G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --partition=chatterjee

/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs04/nilanjan/data/zwang/family/simulation/pgs/data/mergedplink \
--indep-pairwise 1000kb 0.01
