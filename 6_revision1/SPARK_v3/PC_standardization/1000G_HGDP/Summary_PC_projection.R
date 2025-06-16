#PRS standardization
#Family-based study
#Ziqiao Wang
#Feb 27, 2025

library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
#make some PC plots
pc_1KG = fread("~/work/family/data/1000G_HGDP/PCA/1000G_HGDP_projection.sscore",header=T)
library(ggplot2)

colnames(pc_1KG)[5:14] = paste0("PC",1:10)
files = list.files("~/work/data/1000G_HGDP/self_identified_race/")
id=list()
for(i in 1:5){
  id[[i]] = fread(paste0("~/work/data/1000G_HGDP/self_identified_race/",files[i]))
  id[[i]]$Race = gsub("\\..*", "", files[i])
}
id_all = do.call("rbind",id)

pc_1KG_race = merge(pc_1KG,id_all,by.x = "IID",by.y = "V2")
save(pc_1KG_race,file="~/work/family/data/1000G_HGDP/PCA/pc_1KG_race.RData")
p <- ggplot(pc_1KG_race, aes(PC1, PC2))
p <- p + geom_point(aes(colour = factor(Race)))+
scale_color_jco()
ggsave(filename=paste0("PCA_1KG_HGDP.pdf"),
       plot=p, 
       path="~/work/family/data/1000G_HGDP/PCA/",
       width=380, height=300, units="mm", dpi=320)


#combine both dataset together
#make some PC plots for UKB
pc_ukb = fread("~/work/family/data/1000G_HGDP/PCA/UKB/UKB_projection.sscore",header=T)

colnames(pc_ukb)[5:14] = paste0("PC",1:10)

pc_1KG_race$Cohort = "1000G+HGDP"
colnames(pc_1KG_race)[16] = "superclass"
pc_ukb$superclass = "UKB"
pc_ukb$Cohort = "UKB"



pc_1KG_ukb = rbind(pc_1KG_race[,c("IID","PC1","PC2","superclass","Cohort")],pc_ukb[,c("IID","PC1","PC2","superclass","Cohort")])


p1 <- ggplot(pc_1KG_ukb, aes(x = PC1, y = PC2, color = factor(superclass), shape = factor(Cohort), size = factor(Cohort))) + 
  
  geom_point() +
  scale_color_jco() + 
  
  scale_shape_manual(values = c(2,20)) + 
  
  scale_size_manual(values = c(3.5,1)) + theme_bw()

ggsave(filename=paste0("PCA_UKB_1KG.pdf"),
       plot=p1, 
       path="~/work/family/data/1000G_HGDP/PCA/UKB/",
       width=380, height=300, units="mm", dpi=320)



#make some PC plots for SPARK
pc_spark = fread("~/work/family/data/1000G_HGDP/PCA/SPARK_projection_redo.sscore",header=T)

colnames(pc_spark)[5:14] = paste0("PC",1:10)
ancestry = fread("~/work/family/data/autism/v3/pheno/SPARK.iWES_v3.2024_08.ancestry.tsv")
pc_spark_ancestry = merge(pc_spark,ancestry[,c(1,2,3)],by.x = "IID",by.y = "spid")
save(pc_spark_ancestry,file="~/work/family/data/1000G_HGDP/PCA/pc_spark_ancestry_redo.RData")
p <- ggplot(pc_spark_ancestry, aes(PC1, PC2))
p <- p + geom_point(aes(colour = factor(superclass)))+
  scale_color_jco()
ggsave(filename=paste0("PCA_SPARK.pdf"),
       plot=p, 
       path="~/work/family/data/1000G_HGDP/PCA/",
       width=380, height=300, units="mm", dpi=320)


#combine both dataset together
pc_1KG_race$Cohort = "1000G+HGDP"
colnames(pc_1KG_race)[16] = "superclass"
pc_spark_ancestry$Cohort = "SPARK"



pc_1KG_spark = rbind(pc_1KG_race[,c("IID","PC1","PC2","superclass","Cohort")],pc_spark_ancestry[,c("IID","PC1","PC2","superclass","Cohort")])


p1 <- ggplot(pc_1KG_spark, aes(x = PC1, y = PC2, color = factor(superclass), shape = factor(Cohort), size = factor(Cohort))) + 
  
  geom_point() +
  scale_color_jco() + 
  
  scale_shape_manual(values = c(2,20)) + 
  
  scale_size_manual(values = c(3.5,1)) + theme_bw()

ggsave(filename=paste0("PCA_SPARK_1KG_redo.pdf"),
       plot=p1, 
       path="~/work/family/data/1000G_HGDP/PCA/",
       width=380, height=300, units="mm", dpi=320)

## remove mean-effects of ancestry using 1000G+HGDP
prs = data.frame(readRDS("~/work/family/data/1000G_HGDP/PCA/pgs_1KGHGDP_ancestry.rds"))
pc_1KG = fread("~/work/family/data/1000G_HGDP/PCA/1000G_HGDP_projection.sscore",header=T)
colnames(pc_1KG)[5:14] = paste0("PC",1:10)

dat = data.frame(merge(prs,pc_1KG,by = "IID"))

fit_mean =lm(dat[,5] ~ PC1+PC2+PC3+PC4+PC5, data=dat)
pc5_coefs = coef(fit_mean)
pc5_comb = cbind(1, data.matrix(dat)[,paste0('PC', 1:5)]) %*% pc5_coefs
prs_1000G_adj = dat[,5] - as.vector(pc5_comb)
dat$ASD_PGS_adj_mean = prs_1000G_adj  

## remove variance-effects of ancestry
var_1000G_adj = prs_1000G_adj^2
fit_var = lm(var_1000G_adj ~ PC1+PC2+PC3+PC4+PC5, data=dat)
pc5_coefs_var = coef(fit_var)
  
pc5_comb_var = cbind(1, data.matrix(dat)[,paste0('PC', 1:5)]) %*% pc5_coefs_var
prs_1000G_adj_var = prs_1000G_adj / sqrt(pc5_comb_var)
dat$ASD_PGS_adj_sd = prs_1000G_adj_var[,1]
dat=data.frame(dat)
save(dat,pc5_coefs,pc5_coefs_var,file="~/work/family/data/1000G_HGDP/PCA/coef_standardization_prs.RData")

# check calibration by ancestry
library(dplyr)
dat %>%
    group_by(ancestry) %>%
    summarize(raw_mean=mean(`SCORE1_SUM`), raw_sd=sd(`SCORE1_SUM`),
              cal_mean=mean(`ASD_PGS_adj_sd`), cal_sd=sd(`ASD_PGS_adj_sd`),
              cal_mean2=mean(`ASD_PGS_adj_mean`), cal_sd2=sd(`ASD_PGS_adj_mean`))


# ancestry calibration for SPARK
load("~/work/family/data/1000G_HGDP/PCA/coef_standardization_prs.RData")
dat_1KG = dat
load("~/work/family/data/1000G_HGDP/PCA/pc_spark_ancestry_redo.RData")
#load SPARK PRS
prs_spark = readRDS("~/work/family/data/autism/v3/pgs_ASD_final_ancestry.rds")
dat = merge(prs_spark,pc_spark_ancestry,by="IID")
colnames(dat)[7]="superclass"

#match the PRS IDs for the coefficients and SPARK data PRS
#use children super class to define each family
load("~/work/family/data/autism/v3/spark_v3/roles_id_complete_casetrio.RData")
id_clean = c(id_final$subject_sp_id,id_final$biomother_sp_id,id_final$biofather_sp_id)
dat = dat[match(id_clean,dat$IID),]

id_final$ancestry = dat$superclass[match(id_final$subject_sp_id,dat$IID)]
tmp = id_final[,c("subject_sp_id","biomother_sp_id","biofather_sp_id","ancestry")]
dat_superclass_offspring <- melt(setDT(tmp), id.vars = "ancestry", variable.name = "IID")

superclass_offspring = dat_superclass_offspring$ancestry[match(dat$IID,dat_superclass_offspring$value)]



pc5_comb = cbind(1, data.matrix(dat)[,paste0('PC', 1:5)]) %*% pc5_coefs
pc5_comb_var = cbind(1, data.matrix(dat)[,paste0('PC', 1:5)]) %*% pc5_coefs_var

dat$ASD_PGS_adj_sd = as.vector((dat$SCORE1_SUM - pc5_comb) / sqrt(pc5_comb_var))
dat$ASD_PGS_adj_mean=NA

for(i in c("EUR","EAS","AFR","AMR","SAS")){
  
  dat$ASD_PGS_adj_mean[which(superclass_offspring == i)] = as.vector((dat$SCORE1_SUM - pc5_comb))[which(superclass_offspring == i)] / sd(dat_1KG$ASD_PGS_adj_mean[which(dat_1KG$ancestry == i)]) 
  
  
}
summary(dat$ASD_PGS_adj_mean)
sd(dat$ASD_PGS_adj_mean,na.rm = T)


#Standardize using mean/var from 1KG
tmp_adj=rep(NA,dim(dat)[1])
for(i in c("EUR","EAS","AFR","AMR","SAS")){
    
    tmp_adj[which(superclass_offspring == i)] = as.vector((dat$SCORE1_SUM[which(superclass_offspring == i)] - mean(dat_1KG$SCORE1_SUM[which(dat_1KG$ancestry == i)]) )) / sd(dat_1KG$SCORE1_SUM[which(dat_1KG$ancestry == i)]) 
    
}
dat = data.frame(cbind(dat,tmp_adj))
  

colnames(dat)[ncol(dat)] = "ASD_PGS_1KGstand"


save(dat,file="~/work/family/data/autism/v3/PRS_SPARK_PCstandardized.RData")

dat$superclass_offspring = superclass_offspring
# check calibration by ancestry
dat %>%
  group_by(superclass_offspring) %>%
  summarize(raw_mean=mean(`SCORE1_SUM`), raw_sd=sd(`SCORE1_SUM`),
            cal_mean=mean(`ASD_PGS_adj_sd`), cal_sd=sd(`ASD_PGS_adj_sd`),
            cal_mean2=mean(`ASD_PGS_adj_mean`), cal_sd2=sd(`ASD_PGS_adj_mean`),
            cal_mean3=mean(`ASD_PGS_1KGstand`), cal_sd3=sd(`ASD_PGS_1KGstand`))
