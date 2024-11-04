#!/bin/bash
#SBATCH --job-name=family
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:30:00



/dcs04/nilanjan/data/zwang/tools/plink \
--merge-list /dcs04/nilanjan/data/zwang/family/simulation/match_geo/plink/mergelist_ukb.txt \
--freq \
--make-bed --out /dcs04/nilanjan/data/zwang/family/simulation/match_geo/plink/data/mergedplink

