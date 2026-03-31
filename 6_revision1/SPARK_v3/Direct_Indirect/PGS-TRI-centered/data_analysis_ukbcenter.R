#Data Analysis to center the parental indirect effect differences using UKB married male/female data
library(PGS.TRI)
load("./family/data/pleiotropy/redo/pgs_pleiotropy_all_standardized.RData")

#load trio information
load("./family/data/autism/v3/spark_v3/roles_id_complete_casetrio.RData")
id_clean = c(id_final$subject_sp_id,id_final$biomother_sp_id,id_final$biofather_sp_id)
dat_all = dat[match(id_clean,dat$IID),]
dat_all = dat_all[,c(1:52)]
trait=c("edu","schizophrenia","depression","bipolar","adhd","bmi_prive","bipolar1","bipolar2","neuroticism","insomnia","chronotype","asd")

#load UKB mean values
UKB_res_diff = read.csv("./family/data/pleiotropy/UKB/diff_var.csv")
UKB_res_diff_married = read.csv("./family/data/pleiotropy/UKB/diff_var_married.csv")
UKB_res_diff$trait = trait
UKB_res_diff_married$trait = trait

id_final$ancestry = dat$superclass[match(id_final$subject_sp_id,dat$IID)]
id_final$ancestry_m = dat$superclass[match(id_final$biomother_sp_id,dat$IID)]
id_final$ancestry_f = dat$superclass[match(id_final$biofather_sp_id,dat$IID)]
id_final$ancestry_family = paste0(id_final$ancestry,",m=",id_final$ancestry_m,",f=",id_final$ancestry_f)
group = c("EUR","AFR","AMR","EAS","SAS")

#load ASD data
load("./family/data/autism/v3/PRS_SPARK_PCstandardized.RData")
identical(dat$IID,dat_all$IID)
#[1] TRUE
dat = cbind(dat_all,dat[,c(23,24)])

indirect_effect_homo = list()

for(i in 1:length(trait)){
  
  id_final$pgs_child_PC_mean = dat[match(id_final$subject_sp_id,dat$IID),2*i+30]
  id_final$pgs_mother_PC_mean = dat[match(id_final$biomother_sp_id,dat$IID),2*i+30]
  id_final$pgs_father_PC_mean = dat[match(id_final$biofather_sp_id,dat$IID),2*i+30]
  
  
  ###########################################################################
  #repeat the analysis for same ancestry of everyone within each family
  ###########################################################################
  id_same_ancestry = id_final[which(id_final$ancestry == id_final$ancestry_f & id_final$ancestry == id_final$ancestry_m),]

  #Analysis for all population (remove unknown), no standardization
  prs_pgscatalog_complete = id_same_ancestry
  prs_pgscatalog_complete = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry %in% c("AFR", "AMR", "EUR", "EAS", "SAS")),]
  
  ASD_effect=list()
  res_indirect_PC_mean= data.frame(array(0,c(3,4)))
  rownames(res_indirect_PC_mean)= c("EUR (N=12813)","EUR (N=12813), UKB centered (N=337386)","EUR (N=12813), UKB married participants centered (N=250066)")
  
  data.model = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry == "EUR"),]
  ASD_effect[[1]]=PGS.TRI(pgs_offspring = data.model$pgs_child_PC_mean,
                          pgs_mother = data.model$pgs_mother_PC_mean,
                          pgs_father = data.model$pgs_father_PC_mean, 
                          GxE_int = FALSE, smalltriosize = F, parental_indirect = T, parental_diff_ref = UKB_res_diff$diff[i])
  res_indirect_PC_mean[1,] = ASD_effect[[1]]$Coefficients_indirect
  res_indirect_PC_mean[2,] = ASD_effect[[1]]$Coefficients_indirect_corrected
    
  data.model = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry == "EUR"),]
    ASD_effect[[2]]=PGS.TRI(pgs_offspring = data.model$pgs_child_PC_mean,
                            pgs_mother = data.model$pgs_mother_PC_mean,
                            pgs_father = data.model$pgs_father_PC_mean, 
                            GxE_int = FALSE, smalltriosize = F, parental_indirect = T, parental_diff_ref = UKB_res_diff_married$diff[i])
    res_indirect_PC_mean[3,] = ASD_effect[[2]]$Coefficients_indirect_corrected
    
    
  
  colnames(res_indirect_PC_mean)=colnames(ASD_effect[[1]]$Coefficients_indirect)
  res_indirect_PC_mean$standardize= "PC based (with 1KG fixed SD)"
  
  
  indirect_effect_homo[[i]] = res_indirect_PC_mean
  
  
  
}
names(indirect_effect_homo) = trait
save(indirect_effect_homo,file="./family/data/pleiotropy/redo/results_pleiotropy_all_ukbcenter.RData")

library(tidyverse)
indirect_effect_homo_save <-indirect_effect_homo %>% 
  bind_rows(.id = "trait")
colnames(indirect_effect_homo_save)[2:5] = colnames(ASD_effect[[1]]$Coefficients_indirect)
indirect_effect_homo_save$method = rep( c("EUR (N=12813)","EUR (N=12813), UKB centered (N=337386)","EUR (N=12813), UKB married participants centered (N=250066)"),length(trait))
write.csv(indirect_effect_homo_save,file="./family/data/pleiotropy/redo/results_pleiotropy_indirect_homo_ukbcenter.csv",row.names = F)



