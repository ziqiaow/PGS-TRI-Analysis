#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-20:30:00

/dcs04/nilanjan/data/zwang/tools/plink \
--merge-list /dcs04/nilanjan/data/zwang/PRS_autism/plink/data/mergelist.txt --make-bed --out /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/UKB/data/mergedplink

#You can then project onto those PCs with UKB
/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/UKB/data/mergedplink \
--extract /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/data/mergedplink.spark.snp.prune.in \
--read-freq /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/ref_pcs.acount \
--score /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/ref_pcs.eigenvec.allele 2 5 header-read no-mean-imputation cols=+scoresums,-scoreavgs \
--score-col-nums 6-15 \
--out /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/UKB/UKB_projection_nostandardize


#You can then project onto those PCs with UKB, variance standardize
/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/UKB/data/mergedplink \
--extract /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/data/mergedplink.spark.snp.prune.in \
--read-freq /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/ref_pcs.acount \
--score /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/ref_pcs.eigenvec.allele 2 5 header-read no-mean-imputation \
               variance-standardize \
--score-col-nums 6-15 \
--out /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/UKB/UKB_projection
