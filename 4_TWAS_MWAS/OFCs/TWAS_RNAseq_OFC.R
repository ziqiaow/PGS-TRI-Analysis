
#RNAseq
library(data.table)
type="RNAseq"
setwd(paste0("/dcs04/nilanjan/data/zwang/phewas/plink/oc/results/scores/",type))
load("/dcs04/nilanjan/data/zwang/phewas/plink/oc/oc_trio.RData")
source("/users/zwang4/family/R/PGS-TRI.R")

filenames <- list.files(paste0("/dcs04/nilanjan/data/zwang/phewas/plink/oc/results/scores/",type), pattern=".sscore$", full.names=TRUE)
info=fread(paste0("/dcs04/nilanjan/data/zwang/phewas/data/",type,"_trait_validation_results_with_OMICSPRED_ID.csv"))

res_pgs_main=array(0,c(length(filenames),3*6))
res_indirect=array(0,c(length(filenames),3*6))
colnames(res_pgs_main)=colnames(res_indirect)=c("beta_all_all","se_all_all","p_all_all",
                                               "beta_all_EUR","se_all_EUR","p_all_EUR",
                                               "beta_all_AS","se_all_AS","p_all_AS",
                                               "beta_clp_all","se_clp_all","p_clp_all",
                                               "beta_clp_EUR","se_clp_EUR","p_clp_EUR",
                                               "beta_clp_AS","se_clp_AS","p_clp_AS")




for(j in 1:length(filenames)){
  print(j)
  plink2=fread(filenames[j])
  data.model$offspring_prs=plink2$SCORE1_SUM[match(data.model$offspring_id,plink2$IID)]
  data.model$mother_prs=plink2$SCORE1_SUM[match(data.model$mother_id,plink2$IID)]
  data.model$father_prs=plink2$SCORE1_SUM[match(data.model$father_id,plink2$IID)]
  
  
  #standardize the PGS values by race
  data.model$offspring_prs_standardized=0
  data.model$mother_prs_standardized=0
  data.model$father_prs_standardized=0
  
  for(i in c("EUR","ASA")){
    tmp=data.model[which(data.model$proband_race2 == i),]
    prs_all=c(tmp$offspring_prs,tmp$mother_prs,tmp$father_prs)
    tmp$offspring_prs_standardized= (tmp$offspring_prs - mean(prs_all))/sd(prs_all)
    tmp$mother_prs_standardized= (tmp$mother_prs - mean(prs_all))/sd(prs_all)
    tmp$father_prs_standardized= (tmp$father_prs - mean(prs_all))/sd(prs_all)
    
    data.model$offspring_prs_standardized[which(data.model$proband_race2 == i)] = tmp$offspring_prs_standardized
    data.model$mother_prs_standardized[which(data.model$proband_race2 == i)] = tmp$mother_prs_standardized
    data.model$father_prs_standardized[which(data.model$proband_race2 == i)] = tmp$father_prs_standardized
    
    
  }
  
  
  data.model1=data.model[complete.cases(data.model),]
  dim(data.model1)
  
  #all/all
  pgs_main=PGS.TRI(data.model1$offspring_prs_standardized,
                   data.model1$mother_prs_standardized,
                   data.model1$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_indirect = T)
  res_pgs_main[j,1:3]=pgs_main$res_beta[c(1,2,4)]
  res_indirect[j,1:3]=pgs_main$res_delta[c(1,2,4)]
  
  
  
  
  #EUR/all
  data.model2=data.model1[which(data.model1$proband_race2=="EUR" ),]
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_indirect = T)
  res_pgs_main[j,4:6]=pgs_main$res_beta[c(1,2,4)]
  res_indirect[j,4:6]=pgs_main$res_delta[c(1,2,4)]
  
  
  #EUR/all
  data.model2=data.model1[which(data.model1$proband_race2=="ASA" ),]
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_indirect = T)
  res_pgs_main[j,7:9]=pgs_main$res_beta[c(1,2,4)]
  res_indirect[j,7:9]=pgs_main$res_delta[c(1,2,4)]
  
  
  
  
  
  #all/CL/P
  data.model2=data.model1[which(data.model1$clp=="CL/P"),]
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_indirect = T)
  res_pgs_main[j,10:12]=pgs_main$res_beta[c(1,2,4)]
  res_indirect[j,10:12]=pgs_main$res_delta[c(1,2,4)]
  
  
  #EUR, CL/P
  data.model2=data.model1[which(data.model1$proband_race2=="EUR" & data.model1$clp=="CL/P"),]
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_indirect = T)
  res_pgs_main[j,13:15]=pgs_main$res_beta[c(1,2,4)]
  res_indirect[j,13:15]=pgs_main$res_delta[c(1,2,4)]
  
  
  
  
  #Asian, CL/P
  data.model2=data.model1[which(data.model1$proband_race2=="ASA"  & data.model1$clp=="CL/P"),]
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_indirect = T)
  
  res_pgs_main[j,16:18]=pgs_main$res_beta[c(1,2,4)]
  res_indirect[j,16:18]=pgs_main$res_delta[c(1,2,4)]
  
}

name_tmp=gsub(".*/","",filenames)
name_tmp=gsub("\\..*","",name_tmp)
rownames(res_pgs_main)=rownames(res_indirect)=name_tmp
save(res_pgs_main,res_indirect,file=paste0("/dcs04/nilanjan/data/zwang/phewas/plink/oc/results/summary/redo/",type,"_results.RData"))
write.csv(res_pgs_main,file=paste0("/dcs04/nilanjan/data/zwang/phewas/plink/oc/results/summary/redo/",type,"_pgsmain.csv"))
write.csv(res_indirect,file=paste0("/dcs04/nilanjan/data/zwang/phewas/plink/oc/results/summary/redo/",type,"_indirect.csv"))

