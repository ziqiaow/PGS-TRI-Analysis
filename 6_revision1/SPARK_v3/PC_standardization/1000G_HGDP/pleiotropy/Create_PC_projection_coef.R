library(data.table)

trait=c("edu","schizophrenia","depression","bipolar","adhd","bmi_prive","bipolar1","bipolar2","neuroticism","insomnia","chronotype")

#pgs_all=array(0,c(3280,length(trait)))
pgs_all=fread(paste0("/dcs04/nilanjan/data/zwang/family/1000G_HGDP/pleiotropy/edu/score_chr",1,".sscore"),header = T)
pgs_all=pgs_all[,-c(3,4,5)]
h=0
for(i in trait){
  h=h+1
  pgs_each=fread(paste0("/dcs04/nilanjan/data/zwang/family/1000G_HGDP/pleiotropy/",i,"/score_chr",1,".sscore"),header = T)
  score = pgs_each$SCORE1_SUM
  for(j in 2:22){
    tmp=fread(paste0("/dcs04/nilanjan/data/zwang/family/1000G_HGDP/pleiotropy/",i,"/score_chr",j,".sscore"),header = T)
    score_tmp = tmp$SCORE1_SUM
    score = score + score_tmp
  }
  pgs_all = cbind(pgs_all,score)
  
}
colnames(pgs_all)=c("FID","IID",trait)
saveRDS(pgs_all,file="/dcs04/nilanjan/data/zwang/family/1000G_HGDP/pleiotropy/pgs_all_1KGHGDP.rds")

pgs = readRDS("/dcs04/nilanjan/data/zwang/family/1000G_HGDP/pleiotropy/pgs_all_1KGHGDP.rds")
dim(pgs)
#3280 12
 
#load ancestry information
#/dcs05/ladd/NDEpi/data/MasterCohortData/GEARS/SPARK/iWES_v3/metadata/SPARK.iWES_v3.2024_08.ancestry.tsv
setwd("/dcs04/nilanjan/data/ydun/HGDP+1000G/GRCh38/PC_trained_race")
files  <- list.files(pattern = "\\.id$")
ancestry = list()
h=0
for(i in files){
  h=h+1
  ancestry[[h]] = fread(i)
  
}
names(ancestry) = c("AFR","AMR","EAS","EUR","SAS")
pgs$ancestry = 0
pgs$ancestry[match(ancestry[[1]]$V2,pgs$IID)] = "AFR"
pgs$ancestry[match(ancestry[[2]]$V2,pgs$IID)] = "AMR"
pgs$ancestry[match(ancestry[[3]]$V2,pgs$IID)] = "EAS"
pgs$ancestry[match(ancestry[[4]]$V2,pgs$IID)] = "EUR"
pgs$ancestry[match(ancestry[[5]]$V2,pgs$IID)] = "SAS"
saveRDS(pgs,file="/dcs04/nilanjan/data/zwang/family/1000G_HGDP/pleiotropy/pgs_all_1KGHGDP_ancestry.rds")
write.table(pgs,file="/dcs04/nilanjan/data/zwang/family/1000G_HGDP/pleiotropy/pgs_all_1KGHGDP_ancestry.txt",quote = F)





## remove mean-effects of ancestry using 1000G+HGDP
pgs = data.frame(readRDS("/dcs04/nilanjan/data/zwang/family/1000G_HGDP/pleiotropy/pgs_all_1KGHGDP_ancestry.rds"))
pc_1KG = fread("/dcs04/nilanjan/data/zwang/family/1000G_HGDP/plink/PCA/1000G_HGDP_projection.sscore",header=T)
colnames(pc_1KG)[5:14] = paste0("PC",1:10)

dat = data.frame(merge(pgs,pc_1KG,by = "IID"))

coefs_mean = list()
coefs_var = list()
for(i in 1:length(trait)){
  
  fit_mean =lm(dat[,i+2] ~ PC1+PC2+PC3+PC4+PC5, data=dat)
  pc5_coefs = coef(fit_mean)
  pc5_comb = cbind(1, data.matrix(dat)[,paste0('PC', 1:5)]) %*% pc5_coefs
  pgs_1000G_adj = dat[,i+2] - as.vector(pc5_comb)
  
  ## remove variance-effects of ancestry
  var_1000G_adj = pgs_1000G_adj^2
  fit_var = lm(var_1000G_adj ~ PC1+PC2+PC3+PC4+PC5, data=dat)
  pc5_coefs_var = coef(fit_var)
  
  pc5_comb_var = cbind(1, data.matrix(dat)[,paste0('PC', 1:5)]) %*% pc5_coefs_var
  pgs_1000G_adj_var = pgs_1000G_adj / sqrt(pc5_comb_var)
  #rownames(pgs_1000G_adj_var) = dat$IID
  #colnames(pgs_1000G_adj_var) = paste0(trait[i],"_PCstand")
  dat = cbind(dat,pgs_1000G_adj_var[,1],pgs_1000G_adj)
 
  
  coefs_mean[[i]] = pc5_coefs
  coefs_var[[i]] = pc5_coefs_var
}
colnames(dat)[-c(1:27)] = paste0(rep(trait,each=2),c("_adj_sd","_adj_mean"))
names(coefs_mean) = names(coefs_var) = trait
save(dat,coefs_mean,coefs_var,file="/dcs04/nilanjan/data/zwang/family/1000G_HGDP/pleiotropy/coef_standardization_prs.RData")
