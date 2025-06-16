#Pleiotropy traits PGS results summary
#3/3/2025
#first summarize the PRS values
trait=c("edu","schizophrenia","depression","bipolar","adhd","bmi_prive","bipolar1","bipolar2","neuroticism","insomnia","chronotype")
i=1
pgs_all=fread(paste0("/dcs04/nilanjan/data/zwang/family/pleiotropy/redo/",trait[i],"/",trait[i],"_score.sscore"),header = T)
for(i in 2:length(trait)){
  tmp = fread(paste0("/dcs04/nilanjan/data/zwang/family/pleiotropy/redo/",trait[i],"/",trait[i],"_score.sscore"),header = T)
  score_tmp = tmp$SCORE1_SUM
  pgs_all=cbind(pgs_all,score_tmp)
}
colnames(pgs_all)[-c(1:4)] = trait
save(pgs_all,file="/dcs04/nilanjan/data/zwang/family/pleiotropy/redo/pgs_pleiotropy_all.RData")
write.table(pgs_all,file="/dcs04/nilanjan/data/zwang/family/pleiotropy/redo/pgs_pleiotropy_all.txt",quote = F,col.names = T,row.names = F)




#Next standardize using PC from 1KG
# ancestry calibration for SPARK
library(data.table)
load("/dcs04/nilanjan/data/zwang/family/1000G_HGDP/pleiotropy/coef_standardization_prs.RData")
dat_1KG = dat
load("/dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/pc_spark_ancestry_redo.RData")
#load("C:/Users/ziqia/Desktop/work/family/data/1000G_HGDP/PCA/pc_spark_ancestry_redo.RData")
#load SPARK PGSs
load("/dcs04/nilanjan/data/zwang/family/pleiotropy/redo/pgs_pleiotropy_all.RData")
dat = data.frame(merge(pgs_all,pc_spark_ancestry,by="IID"))

#use children super class to define each family
load("/dcs04/nilanjan/data/zwang/family/autism/roles_id_complete_casetrio.RData")
id_clean = c(id_final$subject_sp_id,id_final$biomother_sp_id,id_final$biofather_sp_id)
dat = dat[match(id_clean,dat$IID),]

id_final$ancestry = dat$superclass[match(id_final$subject_sp_id,dat$IID)]
tmp = id_final[,c("subject_sp_id","biomother_sp_id","biofather_sp_id","ancestry")]
dat_superclass_offspring <- melt(setDT(tmp), id.vars = "ancestry", variable.name = "IID")

superclass_offspring = dat_superclass_offspring$ancestry[match(dat$IID,dat_superclass_offspring$value)]

#match the PRS IDs for the coefficients and SPARK data PRS
for(i in 1:length(trait)){
pc5_comb = cbind(1, data.matrix(dat)[,paste0('PC', 1:5)]) %*% coefs_mean[[i]]
pc5_comb_var = cbind(1, data.matrix(dat)[,paste0('PC', 1:5)]) %*% coefs_var[[i]]

tmp_adj = as.vector((dat[,i+4] - pc5_comb) / sqrt(pc5_comb_var))
dat = data.frame(cbind(dat,tmp_adj))

#divided by sd of 1KG adjusted PRS (for each family, based on children's superclass)
tmp_adj=rep(NA,dim(dat)[1])
for(j in c("EUR","EAS","AFR","AMR","SAS")){
  
  tmp_adj[which(superclass_offspring == j)] = as.vector((dat[,i+4] - pc5_comb))[which(superclass_offspring == j)] / sd(dat_1KG[which(dat_1KG$ancestry == j), (2*i+27)]) 
  
}

dat = data.frame(cbind(dat,tmp_adj))

}
colnames(dat)[(ncol(dat)-length(trait)*2+1):ncol(dat)] = paste0(rep(trait,each=2),c("_PCstand_sd","_PCstand_mean"))


#Standardize using mean/var from 1KG
for(i in 1:length(trait)){
  tmp_adj=rep(NA,dim(dat)[1])
for(j in c("EUR","EAS","AFR","AMR","SAS")){
  
  tmp_adj[which(superclass_offspring == j)] = as.vector((dat[which(superclass_offspring == j),i+4] - mean(dat_1KG[which(dat_1KG$ancestry == j), (i+2)]) )) / sd(dat_1KG[which(dat_1KG$ancestry == j), (i+2)]) 
  
}
  dat = data.frame(cbind(dat,tmp_adj))
  
}
colnames(dat)[(ncol(dat)-length(trait)+1):ncol(dat)] = paste0(trait,c("_1KGstand"))

save(dat,file="/dcs04/nilanjan/data/zwang/family/pleiotropy/redo/pgs_pleiotropy_all_standardized.RData")
