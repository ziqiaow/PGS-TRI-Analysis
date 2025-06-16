#Data Analysis
#only use the people in case-parent trios
#load SPARK PRS
#load("/dcs04/nilanjan/data/zwang/family/pleiotropy/redo/pgs_pleiotropy_all_standardized.RData")
load("~/work/family/data/pleiotropy/redo/pgs_pleiotropy_all_standardized.RData")

#load trio information
load("~/work/family/data/autism/v3/spark_v3/roles_id_complete_casetrio.RData")
#load("/dcs04/nilanjan/data/zwang/family/autism/roles_id_complete_casetrio.RData")
id_clean = c(id_final$subject_sp_id,id_final$biomother_sp_id,id_final$biofather_sp_id)
dat = dat[match(id_clean,dat$IID),]
trait=c("edu","schizophrenia","depression","bipolar","adhd","bmi_prive","bipolar1","bipolar2","neuroticism","insomnia","chronotype")


source("~/work/family/R/PGS.TRI/R/PGS-TRI.R")
source("~/work/family/R/PGS.TRI/R/pTDT.R")


id_final$ancestry = dat$superclass[match(id_final$subject_sp_id,dat$IID)]
id_final$ancestry_m = dat$superclass[match(id_final$biomother_sp_id,dat$IID)]
id_final$ancestry_f = dat$superclass[match(id_final$biofather_sp_id,dat$IID)]
id_final$ancestry_family = paste0(id_final$ancestry,",m=",id_final$ancestry_m,",f=",id_final$ancestry_f)


table(id_final$ancestry)
#    AFR      AMR      EAS      EUR      SAS UNKNOWN0 UNKNOWN1 UNKNOWN2 
#1235     2410      442    13668      628      790       19        1 
group = c("EUR","AFR","AMR","EAS","SAS")

direct_effect_all = list()
indirect_effect_all = list()
direct_effect_homo = list()
indirect_effect_homo = list()



for(i in 1:length(trait)){
  
  id_final$pgs_child = dat[match(id_final$subject_sp_id,dat$IID),i+4]
  id_final$pgs_mother = dat[match(id_final$biomother_sp_id,dat$IID),i+4]
  id_final$pgs_father = dat[match(id_final$biofather_sp_id,dat$IID),i+4]
  
  
  id_final$pgs_child_PC_sd = dat[match(id_final$subject_sp_id,dat$IID),2*i+29]
  id_final$pgs_mother_PC_sd = dat[match(id_final$biomother_sp_id,dat$IID),2*i+29]
  id_final$pgs_father_PC_sd = dat[match(id_final$biofather_sp_id,dat$IID),2*i+29]
  
  id_final$pgs_child_PC_mean = dat[match(id_final$subject_sp_id,dat$IID),2*i+30]
  id_final$pgs_mother_PC_mean = dat[match(id_final$biomother_sp_id,dat$IID),2*i+30]
  id_final$pgs_father_PC_mean = dat[match(id_final$biofather_sp_id,dat$IID),2*i+30]
  
  
  id_final$pgs_child_1KG = dat[match(id_final$subject_sp_id,dat$IID),i+52]
  id_final$pgs_mother_1KG = dat[match(id_final$biomother_sp_id,dat$IID),i+52]
  id_final$pgs_father_1KG = dat[match(id_final$biofather_sp_id,dat$IID),i+52]
  
  
  #Analysis for all population (including unknown), no standardization
  prs_pgscatalog_complete = id_final
  ASD_effect=list()
  res_direct = res_indirect= data.frame(array(0,c(6,4)))
  rownames(res_direct) =rownames(res_indirect)= c("All (N=19193)","EUR (N=13668)","AFR (N=1235)","AMR (N=2410)",
                                                        "EAS (N=442)","SAS (N=628)")
  
  ASD_effect[[1]]=PGS.TRI(pgs_offspring = prs_pgscatalog_complete$pgs_child,
                          pgs_mother = prs_pgscatalog_complete$pgs_mother,
                          pgs_father = prs_pgscatalog_complete$pgs_father, 
                          GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
  res_direct[1,] = ASD_effect[[1]]$res_beta
  res_indirect[1,] = ASD_effect[[1]]$res_delta
  
  for(j in 2:6){
    
    data.model = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry == group[j-1]),]
    ASD_effect[[j]]=PGS.TRI(pgs_offspring = data.model$pgs_child,
                            pgs_mother = data.model$pgs_mother,
                            pgs_father = data.model$pgs_father, 
                            GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
    res_direct[j,] = ASD_effect[[j]]$res_beta
    res_indirect[j,] = ASD_effect[[j]]$res_delta
    
    
  }
  
  colnames(res_direct)=colnames(res_indirect)=colnames(ASD_effect[[1]]$res_beta)
  res_direct$standardize =res_indirect$standardize= "No standardization"
  
  
  #Analysis for PC standardized, all population (including unknown)
  ASD_effect=list()
  res_direct_PC = res_indirect_PC= data.frame(array(0,c(6,4)))
  rownames(res_direct_PC) =rownames(res_indirect_PC)= c("All (N=19193)","EUR (N=13668)","AFR (N=1235)","AMR (N=2410)",
                                                        "EAS (N=442)","SAS (N=628)")
  
  ASD_effect[[1]]=PGS.TRI(pgs_offspring = prs_pgscatalog_complete$pgs_child_PC_sd,
                          pgs_mother = prs_pgscatalog_complete$pgs_mother_PC_sd,
                          pgs_father = prs_pgscatalog_complete$pgs_father_PC_sd, 
                          GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
  res_direct_PC[1,] = ASD_effect[[1]]$res_beta
  res_indirect_PC[1,] = ASD_effect[[1]]$res_delta
  
  for(j in 2:6){
    
    data.model = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry == group[j-1]),]
    ASD_effect[[j]]=PGS.TRI(pgs_offspring = data.model$pgs_child_PC_sd,
                            pgs_mother = data.model$pgs_mother_PC_sd,
                            pgs_father = data.model$pgs_father_PC_sd, 
                            GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
    res_direct_PC[j,] = ASD_effect[[j]]$res_beta
    res_indirect_PC[j,] = ASD_effect[[j]]$res_delta
    
    
  }

  colnames(res_direct_PC)=colnames(res_indirect_PC)=colnames(ASD_effect[[1]]$res_beta)
  res_direct_PC$standardize =res_indirect_PC$standardize= "PC based (with SD)"
  
  
  #Analysis for PC mean standardized, excluding unknown
  prs_pgscatalog_complete = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry %in% c("AFR", "AMR", "EUR", "EAS", "SAS")),]
  
  ASD_effect=list()
  res_direct_PC_mean = res_indirect_PC_mean= data.frame(array(0,c(6,4)))
  rownames(res_direct_PC_mean) =rownames(res_indirect_PC_mean)= c("All (N=18383)","EUR (N=13668)","AFR (N=1235)","AMR (N=2410)",
                                                                "EAS (N=442)","SAS (N=628)")
  
  ASD_effect[[1]]=PGS.TRI(pgs_offspring = prs_pgscatalog_complete$pgs_child_PC_mean,
                          pgs_mother = prs_pgscatalog_complete$pgs_mother_PC_mean,
                          pgs_father = prs_pgscatalog_complete$pgs_father_PC_mean, 
                          GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
  res_direct_PC_mean[1,] = ASD_effect[[1]]$res_beta
  res_indirect_PC_mean[1,] = ASD_effect[[1]]$res_delta
  
  for(j in 2:6){
    
    data.model = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry == group[j-1]),]
    ASD_effect[[j]]=PGS.TRI(pgs_offspring = data.model$pgs_child_PC_mean,
                            pgs_mother = data.model$pgs_mother_PC_mean,
                            pgs_father = data.model$pgs_father_PC_mean, 
                            GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
    res_direct_PC_mean[j,] = ASD_effect[[j]]$res_beta
    res_indirect_PC_mean[j,] = ASD_effect[[j]]$res_delta
    
    
  }
  
  colnames(res_direct_PC_mean)=colnames(res_indirect_PC_mean)=colnames(ASD_effect[[1]]$res_beta)
  res_direct_PC_mean$standardize =res_indirect_PC_mean$standardize= "PC based (with 1KG fixed SD)"
  
  
  #Analysis for 1KG mean/sd standardized, excluding unknown
  ASD_effect=list()
  res_direct_1KG = res_indirect_1KG= data.frame(array(0,c(6,4)))
  rownames(res_direct_1KG) =rownames(res_indirect_1KG)= c("All (N=18383)","EUR (N=13668)","AFR (N=1235)","AMR (N=2410)",
                                                                  "EAS (N=442)","SAS (N=628)")
  
  ASD_effect[[1]]=PGS.TRI(pgs_offspring = prs_pgscatalog_complete$pgs_child_1KG,
                          pgs_mother = prs_pgscatalog_complete$pgs_mother_1KG,
                          pgs_father = prs_pgscatalog_complete$pgs_father_1KG, 
                          GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
  res_direct_1KG[1,] = ASD_effect[[1]]$res_beta
  res_indirect_1KG[1,] = ASD_effect[[1]]$res_delta
  
  for(j in 2:6){
    
    data.model = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry == group[j-1]),]
    ASD_effect[[j]]=PGS.TRI(pgs_offspring = data.model$pgs_child_1KG,
                            pgs_mother = data.model$pgs_mother_1KG,
                            pgs_father = data.model$pgs_father_1KG, 
                            GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
    res_direct_1KG[j,] = ASD_effect[[j]]$res_beta
    res_indirect_1KG[j,] = ASD_effect[[j]]$res_delta
    
    
  }
  
  colnames(res_direct_1KG)=colnames(res_indirect_1KG)=colnames(ASD_effect[[1]]$res_beta)
  res_direct_1KG$standardize =res_indirect_1KG$standardize= "1KG Fixed Mean/SD based"
  
  
  direct_effect_all[[i]] = rbind(res_direct,res_direct_PC,res_direct_PC_mean,res_direct_1KG)
  indirect_effect_all[[i]] = rbind(res_indirect,res_indirect_PC,res_indirect_PC_mean,res_indirect_1KG)
 
  
  
  
  
  ###########################################################################
  #repeat the analysis for same ancestry of everyone within each family
  ###########################################################################
  id_same_ancestry = id_final[which(id_final$ancestry == id_final$ancestry_f & id_final$ancestry == id_final$ancestry_m),]
  
  #Analysis for all population (remove unknown), no standardization
  prs_pgscatalog_complete = id_same_ancestry
  prs_pgscatalog_complete = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry %in% c("AFR", "AMR", "EUR", "EAS", "SAS")),]
  
  ASD_effect=list()
  res_direct = res_indirect= data.frame(array(0,c(6,4)))
  rownames(res_direct) =rownames(res_indirect)= c("All (N=15876)","EUR (N=12813)","AFR (N=792)","AMR (N=1302)",
                                                  "EAS (N=415)","SAS (N=554)")
  
  ASD_effect[[1]]=PGS.TRI(pgs_offspring = prs_pgscatalog_complete$pgs_child,
                          pgs_mother = prs_pgscatalog_complete$pgs_mother,
                          pgs_father = prs_pgscatalog_complete$pgs_father, 
                          GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
  res_direct[1,] = ASD_effect[[1]]$res_beta
  res_indirect[1,] = ASD_effect[[1]]$res_delta
  
  for(j in 2:6){
    
    data.model = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry == group[j-1]),]
    ASD_effect[[j]]=PGS.TRI(pgs_offspring = data.model$pgs_child,
                            pgs_mother = data.model$pgs_mother,
                            pgs_father = data.model$pgs_father, 
                            GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
    res_direct[j,] = ASD_effect[[j]]$res_beta
    res_indirect[j,] = ASD_effect[[j]]$res_delta
    
    
  }
  
  colnames(res_direct)=colnames(res_indirect)=colnames(ASD_effect[[1]]$res_beta)
  res_direct$standardize =res_indirect$standardize= "No standardization"
  
  
  #Analysis for PC standardized, all population (remove unknown)
  ASD_effect=list()
  res_direct_PC = res_indirect_PC= data.frame(array(0,c(6,4)))
  rownames(res_direct_PC) =rownames(res_indirect_PC)= c("All (N=15876)","EUR (N=12813)","AFR (N=792)","AMR (N=1302)",
                                                        "EAS (N=415)","SAS (N=554)")
  
  ASD_effect[[1]]=PGS.TRI(pgs_offspring = prs_pgscatalog_complete$pgs_child_PC_sd,
                          pgs_mother = prs_pgscatalog_complete$pgs_mother_PC_sd,
                          pgs_father = prs_pgscatalog_complete$pgs_father_PC_sd, 
                          GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
  res_direct_PC[1,] = ASD_effect[[1]]$res_beta
  res_indirect_PC[1,] = ASD_effect[[1]]$res_delta
  
  for(j in 2:6){
    
    data.model = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry == group[j-1]),]
    ASD_effect[[j]]=PGS.TRI(pgs_offspring = data.model$pgs_child_PC_sd,
                            pgs_mother = data.model$pgs_mother_PC_sd,
                            pgs_father = data.model$pgs_father_PC_sd, 
                            GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
    res_direct_PC[j,] = ASD_effect[[j]]$res_beta
    res_indirect_PC[j,] = ASD_effect[[j]]$res_delta
    
    
  }
  
  colnames(res_direct_PC)=colnames(res_indirect_PC)=colnames(ASD_effect[[1]]$res_beta)
  res_direct_PC$standardize =res_indirect_PC$standardize= "PC based (with SD)"
  
  
  #Analysis for PC mean standardized, excluding unknown
  prs_pgscatalog_complete = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry %in% c("AFR", "AMR", "EUR", "EAS", "SAS")),]
  
  ASD_effect=list()
  res_direct_PC_mean = res_indirect_PC_mean= data.frame(array(0,c(6,4)))
  rownames(res_direct_PC_mean) =rownames(res_indirect_PC_mean)= c("All (N=15876)","EUR (N=12813)","AFR (N=792)","AMR (N=1302)",
                                                                  "EAS (N=415)","SAS (N=554)")
  
  ASD_effect[[1]]=PGS.TRI(pgs_offspring = prs_pgscatalog_complete$pgs_child_PC_mean,
                          pgs_mother = prs_pgscatalog_complete$pgs_mother_PC_mean,
                          pgs_father = prs_pgscatalog_complete$pgs_father_PC_mean, 
                          GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
  res_direct_PC_mean[1,] = ASD_effect[[1]]$res_beta
  res_indirect_PC_mean[1,] = ASD_effect[[1]]$res_delta
  
  for(j in 2:6){
    
    data.model = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry == group[j-1]),]
    ASD_effect[[j]]=PGS.TRI(pgs_offspring = data.model$pgs_child_PC_mean,
                            pgs_mother = data.model$pgs_mother_PC_mean,
                            pgs_father = data.model$pgs_father_PC_mean, 
                            GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
    res_direct_PC_mean[j,] = ASD_effect[[j]]$res_beta
    res_indirect_PC_mean[j,] = ASD_effect[[j]]$res_delta
    
    
  }
  
  colnames(res_direct_PC_mean)=colnames(res_indirect_PC_mean)=colnames(ASD_effect[[1]]$res_beta)
  res_direct_PC_mean$standardize =res_indirect_PC_mean$standardize= "PC based (with 1KG fixed SD)"
  
  
  #Analysis for 1KG mean/sd standardized, excluding unknown
  ASD_effect=list()
  res_direct_1KG = res_indirect_1KG= data.frame(array(0,c(6,4)))
  rownames(res_direct_1KG) =rownames(res_indirect_1KG)= c("All (N=15876)","EUR (N=12813)","AFR (N=792)","AMR (N=1302)",
                                                          "EAS (N=415)","SAS (N=554)")
  
  ASD_effect[[1]]=PGS.TRI(pgs_offspring = prs_pgscatalog_complete$pgs_child_1KG,
                          pgs_mother = prs_pgscatalog_complete$pgs_mother_1KG,
                          pgs_father = prs_pgscatalog_complete$pgs_father_1KG, 
                          GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
  res_direct_1KG[1,] = ASD_effect[[1]]$res_beta
  res_indirect_1KG[1,] = ASD_effect[[1]]$res_delta
  
  for(j in 2:6){
    
    data.model = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry == group[j-1]),]
    ASD_effect[[j]]=PGS.TRI(pgs_offspring = data.model$pgs_child_1KG,
                            pgs_mother = data.model$pgs_mother_1KG,
                            pgs_father = data.model$pgs_father_1KG, 
                            GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
    res_direct_1KG[j,] = ASD_effect[[j]]$res_beta
    res_indirect_1KG[j,] = ASD_effect[[j]]$res_delta
    
    
  }
  
  colnames(res_direct_1KG)=colnames(res_indirect_1KG)=colnames(ASD_effect[[1]]$res_beta)
  res_direct_1KG$standardize =res_indirect_1KG$standardize= "1KG Fixed Mean/SD based"
  
  
  direct_effect_homo[[i]] = rbind(res_direct,res_direct_PC,res_direct_PC_mean,res_direct_1KG)
  indirect_effect_homo[[i]] = rbind(res_indirect,res_indirect_PC,res_indirect_PC_mean,res_indirect_1KG)
  
  
  
}
names(direct_effect_all) = names(indirect_effect_all) = names(direct_effect_homo) = names(indirect_effect_homo) = trait
save(direct_effect_all,direct_effect_homo,indirect_effect_all,indirect_effect_homo,file="~/work/family/data/pleiotropy/redo/results_pleiotropy_all.RData")

library(tidyverse)
direct_effect_all_save <-direct_effect_all %>% 
  bind_rows(.id = "trait")

direct_effect_homo_save <-direct_effect_homo %>% 
  bind_rows(.id = "trait")
indirect_effect_all_save <-indirect_effect_all %>% 
  bind_rows(.id = "trait")
indirect_effect_homo_save <-indirect_effect_homo %>% 
  bind_rows(.id = "trait")

write.csv(direct_effect_all_save,file="~/work/family/data/pleiotropy/redo/results_pleiotropy_direct_all.csv")
write.csv(direct_effect_homo_save,file="~/work/family/data/pleiotropy/redo/results_pleiotropy_direct_homo.csv")
write.csv(indirect_effect_all_save,file="~/work/family/data/pleiotropy/redo/results_pleiotropy_indirect_all.csv")
write.csv(indirect_effect_homo_save,file="~/work/family/data/pleiotropy/redo/results_pleiotropy_indirect_homo.csv")


