#!/bin/bash
#SBATCH --job-name=bc
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6



module load plink/1.90b
plink --bfile /dcs04/nilanjan/data/zwang/family/simulation/match_geo/plink/data/mergedplink \
--extract /dcs04/nilanjan/data/zwang/family/simulation/pgs/ea_score_parents_pruned.sscore.vars \
--recode A --out /dcs04/nilanjan/data/zwang/family/simulation/match_geo/plink/data/matched_parents_pruned

