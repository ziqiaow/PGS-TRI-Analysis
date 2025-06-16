#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --time=0-20:30:00


#Combine 1000G+HGDP by chromosomes
/dcs04/nilanjan/data/zwang/tools/plink \
--merge-list /dcs04/nilanjan/data/zwang/PRS_autism/plink/data/mergelist.txt --make-bed --out /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/data/mergedplink



#find SPARK and 1KG overlapping snps and prune for independent SNPs to prepare for PCA calculation
/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/data/mergedplink \
--extract /dcs05/ladd/NDEpi/data/Projects/InProgress/ziqiao/data/mergefinal.bim \
--make-bed --out /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/data/mergedplink.spark.snp \
--indep-pairwise 5000kb 0.1 \

#Calculate PCs
/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/data/mergedplink.spark.snp \
--extract /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/data/mergedplink.spark.snp.prune.in \
--freq counts \
--pca allele-wts vcols=chrom,ref,alt \
--out /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/ref_pcs

#Project onto those PCs with 1000G+HGDP itself
/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/data/mergedplink.spark.snp \
--extract /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/data/mergedplink.spark.snp.prune.in \
--read-freq /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/ref_pcs.acount \
--score /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/ref_pcs.eigenvec.allele 2 5 header-read no-mean-imputation \
               variance-standardize \
--score-col-nums 6-15 \
--out /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/1000G_HGDP_projection



#You can then project onto those PCs with SPARK data
/dcs04/nilanjan/data/zwang/tools/20240105/plink2 \
--bfile /dcs05/ladd/NDEpi/data/Projects/InProgress/ziqiao/data/mergefinal \
--extract /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/data/mergedplink.spark.snp.prune.in \
--read-freq /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/ref_pcs.acount \
--score /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/ref_pcs.eigenvec.allele 2 5 header-read no-mean-imputation \
               variance-standardize \
--score-col-nums 6-15 \
--out /dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/SPARK_projection_redo
