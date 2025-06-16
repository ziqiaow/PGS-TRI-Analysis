#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22
#SBATCH --time=0-09:30:00

#You can then project onto those PCs with UKB data
/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--pfile /dcl01/arking/data/UK_Biobank/static/Genetic/GenotypingArray/Imputed/plink/chr$SLURM_ARRAY_TASK_ID \
--rm-dup exclude-all \
--extract /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/data/mergedplink.spark.snp.prune.in \
--make-bed --out /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/UKB/data/chr$SLURM_ARRAY_TASK_ID

