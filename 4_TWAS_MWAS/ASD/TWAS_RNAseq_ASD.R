#RNAseq
library(data.table)
type="RNAseq"
setwd(paste0("/dcs04/nilanjan/data/zwang/phewas/plink/results/scores/",type))
load("/dcs04/nilanjan/data/zwang/phewas/plink/results/summary_asd_trio.RData")
source("/users/zwang4/family/R/PGS-TRI.R")

filenames <- list.files(paste0("/dcs04/nilanjan/data/zwang/phewas/plink/results/scores/",type), pattern=".sscore$", full.names=TRUE)

info=fread(paste0("/dcs04/nilanjan/data/zwang/phewas/data/",type,"_trait_validation_results_with_OMICSPRED_ID.csv"))

res_pgs_main=array(0,c(length(filenames),3*6))
res_nurture=array(0,c(length(filenames),3*6))
colnames(res_pgs_main)=colnames(res_nurture)=c("beta_all","se_all","p_all",
                         "beta_all_noAMR","se_all_noAMR","p_all_noAMR",
                         "beta_EUR","se_EUR","p_EUR",
                         "beta_AFR","se_AFR","p_AFR",
                         "beta_AMR","se_AMR","p_AMR",
                         "beta_AS","se_AS","p_AS")




for(j in 1:length(filenames)){
  print(j)
plink2=fread(filenames[j])
prs_pgscatalog_complete$offspring_prs=plink2$SCORE1_SUM[match(prs_pgscatalog_complete$offspring_id,plink2$IID)]
prs_pgscatalog_complete$mother_prs=plink2$SCORE1_SUM[match(prs_pgscatalog_complete$mother_id,plink2$IID)]
prs_pgscatalog_complete$father_prs=plink2$SCORE1_SUM[match(prs_pgscatalog_complete$father_id,plink2$IID)]


#standardize the PGS values by race
prs_pgscatalog_complete$offspring_prs_standardized=0
prs_pgscatalog_complete$mother_prs_standardized=0
prs_pgscatalog_complete$father_prs_standardized=0

for(i in unique(prs_pgscatalog_complete$race)){
  tmp=prs_pgscatalog_complete[which(prs_pgscatalog_complete$race == i),]
  prs_all=c(tmp$offspring_prs,tmp$mother_prs,tmp$father_prs)
  tmp$offspring_prs_standardized= (tmp$offspring_prs - mean(prs_all))/sd(prs_all)
  tmp$mother_prs_standardized= (tmp$mother_prs - mean(prs_all))/sd(prs_all)
  tmp$father_prs_standardized= (tmp$father_prs - mean(prs_all))/sd(prs_all)
  
  prs_pgscatalog_complete$offspring_prs_standardized[which(prs_pgscatalog_complete$race == i)] = tmp$offspring_prs_standardized
  prs_pgscatalog_complete$mother_prs_standardized[which(prs_pgscatalog_complete$race == i)] = tmp$mother_prs_standardized
  prs_pgscatalog_complete$father_prs_standardized[which(prs_pgscatalog_complete$race == i)] = tmp$father_prs_standardized
  
  
}


data.model1=prs_pgscatalog_complete[complete.cases(prs_pgscatalog_complete),]
dim(data.model1)

#all
pgs_main=PGS.TRI(data.model1$offspring_prs_standardized,
                 data.model1$mother_prs_standardized,
                 data.model1$father_prs_standardized, 
                 GxE_int = FALSE,
                 parental_indirect = T)
res_pgs_main[j,1:3]=pgs_main$res_beta[c(1,2,4)]
res_nurture[j,1:3]=pgs_main$res_delta[c(1,2,4)]

#all no AMR
data.model2=data.model1[-which(data.model1$race=="AMR"),]
pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                 data.model2$mother_prs_standardized,
                 data.model2$father_prs_standardized, 
                 GxE_int = FALSE,
                 parental_indirect = T)
res_pgs_main[j,4:6]=pgs_main$res_beta[c(1,2,4)]
res_nurture[j,4:6]=pgs_main$res_delta[c(1,2,4)]

#EUR
data.model2=data.model1[which(data.model1$race=="EUR"),]
pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                 data.model2$mother_prs_standardized,
                 data.model2$father_prs_standardized, 
                 GxE_int = FALSE,
                 parental_indirect = T)
res_pgs_main[j,7:9]=pgs_main$res_beta[c(1,2,4)]
res_nurture[j,7:9]=pgs_main$res_delta[c(1,2,4)]


#AFR
data.model2=data.model1[which(data.model1$race=="AFR"),]
pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                 data.model2$mother_prs_standardized,
                 data.model2$father_prs_standardized, 
                 GxE_int = FALSE,
                 parental_indirect = T, smalltriosize = T)
res_pgs_main[j,10:12]=pgs_main$res_beta[c(1,2,4)]
res_nurture[j,10:12]=pgs_main$res_delta[c(1,2,4)]


#AMR
data.model2=data.model1[which(data.model1$race=="AMR"),]
pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                 data.model2$mother_prs_standardized,
                 data.model2$father_prs_standardized, 
                 GxE_int = FALSE,
                 parental_indirect = T)
res_pgs_main[j,13:15]=pgs_main$res_beta[c(1,2,4)]
res_nurture[j,13:15]=pgs_main$res_delta[c(1,2,4)]



#EAS/SAS
data.model2=data.model1[which(data.model1$race=="EAS" | data.model1$race=="SAS"),]
pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                 data.model2$mother_prs_standardized,
                 data.model2$father_prs_standardized, 
                 GxE_int = FALSE,
                 parental_indirect = T, smalltriosize = T)

res_pgs_main[j,16:18]=pgs_main$res_beta[c(1,2,4)]
res_nurture[j,16:18]=pgs_main$res_delta[c(1,2,4)]

}

name_tmp=gsub(".*/","",filenames)
name_tmp=gsub("\\..*","",name_tmp)
rownames(res_pgs_main)=rownames(res_nurture)=name_tmp
save(res_pgs_main,res_nurture,file=paste0("/dcs04/nilanjan/data/zwang/phewas/plink/results/summary/redo/",type,"_results_smallsample.RData"))
write.csv(res_pgs_main,file=paste0("/dcs04/nilanjan/data/zwang/phewas/plink/results/summary/redo/",type,"_pgsmain_smallsample.csv"))
write.csv(res_nurture,file=paste0("/dcs04/nilanjan/data/zwang/phewas/plink/results/summary/redo/",type,"_nurture_smallsample.csv"))