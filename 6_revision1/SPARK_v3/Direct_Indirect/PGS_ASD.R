
library(data.table)
pgs_each=fread(paste0("/dcs04/nilanjan/data/zwang/family/autism/plink/pgs/redo/ASD_score_final.sscore"),header = T)

#load ancestry information
#/dcs05/ladd/NDEpi/data/MasterCohortData/GEARS/SPARK/iWES_v3/metadata/SPARK.iWES_v3.2024_08.ancestry.tsv
ancestry = fread("/dcs05/ladd/NDEpi/data/MasterCohortData/GEARS/SPARK/iWES_v3/metadata/SPARK.iWES_v3.2024_08.ancestry.tsv")
pgs_ancestry = merge(pgs_each,ancestry[,c(1,2,3)],by.x="IID",by.y="spid")
dim(pgs_ancestry)
#[1] 141368      7

saveRDS(pgs_ancestry,file="/dcs04/nilanjan/data/zwang/family/autism/plink/pgs/redo/pgs_ASD_final_ancestry.rds")



source("C:/Users/ziqia/Desktop/work/family/R/PGS.TRI/R/PGS-TRI.R")
source("C:/Users/ziqia/Desktop/work/family/R/PGS.TRI/R/pTDT.R")

#PC normalization results is available in
#this step was done in C:/Users/ziqia/Desktop/work/family/data/1000G_HGDP/PCA/summary_pca.R
load("C:/Users/ziqia/Desktop/work/family/data/autism/v3/PRS_SPARK_PCstandardized.RData")


#load ancestry information
#/dcs05/ladd/NDEpi/data/MasterCohortData/GEARS/SPARK/iWES_v3/metadata/SPARK.iWES_v3.2024_08.ancestry.tsv
#ancestry = fread("/dcs05/ladd/NDEpi/data/MasterCohortData/GEARS/SPARK/iWES_v3/metadata/SPARK.iWES_v3.2024_08.ancestry.tsv")
#ancestry = fread("C:/Users/ziqia/Desktop/work/family/data/autism/v3/pheno/SPARK.iWES_v3.2024_08.ancestry.tsv")
#sample_info = fread("C:/Users/ziqia/Desktop/work/family/data/autism/v3/pheno/SPARK.iWES_v3.2024_08.sample_metadata.tsv")
#asd_info = sample_info[which(sample_info$asd == "TRUE"),]
#id = intersect(asd_info$spid, id_final$subject_sp_id) #19193, this is consistent

#id = intersect(id_final$subject_sp_id,dat$IID)

# asd_trio = asd_info[match(id_final$subject_sp_id,asd_info$spid),]
# identical(asd_trio$father,id_final$biofather_sp_id) #[1] TRUE
# identical(asd_trio$mother,id_final$biomother_sp_id) #[1] TRUE

load("C:/Users/ziqia/Desktop/work/family/data/autism/v3/spark_v3/roles_id_complete_casetrio.RData")


id_final$ancestry = dat$superclass[match(id_final$subject_sp_id,dat$IID)]
id_final$ancestry_m = dat$superclass[match(id_final$biomother_sp_id,dat$IID)]
id_final$ancestry_f = dat$superclass[match(id_final$biofather_sp_id,dat$IID)]
id_final$ancestry_family = paste0(id_final$ancestry,",m=",id_final$ancestry_m,",f=",id_final$ancestry_f)


table(id_final$ancestry)
#    AFR      AMR      EAS      EUR      SAS UNKNOWN0 UNKNOWN1 UNKNOWN2 
#1235     2410      442    13668      628      790       19        1 
group = c("EUR","AFR","AMR","EAS","SAS")


  
  id_final$pgs_child = dat$SCORE1_SUM[match(id_final$subject_sp_id,dat$IID)]
  id_final$pgs_mother = dat$SCORE1_SUM[match(id_final$biomother_sp_id,dat$IID)]
  id_final$pgs_father = dat$SCORE1_SUM[match(id_final$biofather_sp_id,dat$IID)]
  
  
  id_final$pgs_child_PC_sd = dat$ASD_PGS_adj_sd[match(id_final$subject_sp_id,dat$IID)]
  id_final$pgs_mother_PC_sd = dat$ASD_PGS_adj_sd[match(id_final$biomother_sp_id,dat$IID)]
  id_final$pgs_father_PC_sd = dat$ASD_PGS_adj_sd[match(id_final$biofather_sp_id,dat$IID)]
  
  id_final$pgs_child_PC_mean = dat$ASD_PGS_adj_mean[match(id_final$subject_sp_id,dat$IID)]
  id_final$pgs_mother_PC_mean = dat$ASD_PGS_adj_mean[match(id_final$biomother_sp_id,dat$IID)]
  id_final$pgs_father_PC_mean = dat$ASD_PGS_adj_mean[match(id_final$biofather_sp_id,dat$IID)]
  
  
  id_final$pgs_child_1KG = dat$ASD_PGS_1KGstand[match(id_final$subject_sp_id,dat$IID)]
  id_final$pgs_mother_1KG = dat$ASD_PGS_1KGstand[match(id_final$biomother_sp_id,dat$IID)]
  id_final$pgs_father_1KG = dat$ASD_PGS_1KGstand[match(id_final$biofather_sp_id,dat$IID)]
  
  
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
  
  
  direct_effect_all = rbind(res_direct,res_direct_PC,res_direct_PC_mean,res_direct_1KG)
  indirect_effect_all = rbind(res_indirect,res_indirect_PC,res_indirect_PC_mean,res_indirect_1KG)
  
  
  
  
  
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
  write.csv(res_direct_PC_mean,file="C:/Users/ziqia/Desktop/work/family/data/autism/v3/direct_results.csv")
  write.csv(res_indirect_PC_mean,file="C:/Users/ziqia/Desktop/work/family/data/autism/v3/indirect_results.csv")
  
  
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
  
  
direct_effect_homo = rbind(res_direct,res_direct_PC,res_direct_PC_mean,res_direct_1KG)
indirect_effect_homo= rbind(res_indirect,res_indirect_PC,res_indirect_PC_mean,res_indirect_1KG)
  


#pTDT test
res_ptdt = data.frame(array(0,c(6,4)))
rownames(res_ptdt) = c("All (N=15876)","EUR (N=12813)","AFR (N=792)","AMR (N=1302)",
                                                        "EAS (N=415)","SAS (N=554)")

ptdt.test=list()
ptdt.test[[1]]=ptdt(prs_pgscatalog_complete$pgs_child_PC_mean,
                    prs_pgscatalog_complete$pgs_mother_PC_mean,
                    prs_pgscatalog_complete$pgs_father_PC_mean)
res_ptdt[1,]=ptdt.test[[1]]$res_beta

for(j in 2:6){
  
  data.model = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry == group[j-1]),]
  ptdt.test[[j]]=ptdt(data.model$pgs_child_PC_mean,
                      data.model$pgs_mother_PC_mean,
                      data.model$pgs_father_PC_mean)
  res_ptdt[j,]=ptdt.test[[j]]$res_beta
  
}

colnames(res_ptdt)=colnames(ptdt.test[[1]]$res_beta)
write.csv(res_ptdt,file="C:/Users/ziqia/Desktop/work/family/data/autism/v3/ptdt_results.csv")


###########################################################################
#repeat the analysis for different ancestry of everyone within each family
###########################################################################
id_heter_ancestry = id_final[which(id_final$ancestry != id_final$ancestry_f | id_final$ancestry != id_final$ancestry_m | id_final$ancestry_f != id_final$ancestry_m),]

#Analysis for all population (remove unknown), no standardization
prs_pgscatalog_complete = id_heter_ancestry
table(prs_pgscatalog_complete$ancestry)

#  AFR      AMR      EAS      EUR      SAS UNKNOWN0 UNKNOWN1 UNKNOWN2 
#443     1108       27      855       74      765       17        1 



ASD_effect=list()
res_direct = res_indirect= data.frame(array(0,c(7,4)))
rownames(res_direct) =rownames(res_indirect)=  c("All (N=3290)","EUR (N=855)","AFR (N=443)","AMR (N=1108)",
                                                 "EAS (N=27)","SAS (N=74)","All, no unknown (N=2507)")

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
data.model = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry %in% c("AFR", "AMR", "EUR", "EAS", "SAS")),]

ASD_effect[[7]]=PGS.TRI(pgs_offspring = data.model$pgs_child,
                        pgs_mother = data.model$pgs_mother,
                        pgs_father = data.model$pgs_father, 
                        GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
res_direct[7,] = ASD_effect[[7]]$res_beta
res_indirect[7,] = ASD_effect[[7]]$res_delta


colnames(res_direct)=colnames(res_indirect)=colnames(ASD_effect[[1]]$res_beta)
res_direct$standardize =res_indirect$standardize= "No standardization"


#Analysis for PC standardized, all population (include unknown)
ASD_effect=list()
res_direct_PC = res_indirect_PC= data.frame(array(0,c(7,4)))
rownames(res_direct_PC) =rownames(res_indirect_PC)= c("All (N=3290)","EUR (N=855)","AFR (N=443)","AMR (N=1108)",
                                                      "EAS (N=27)","SAS (N=74)","All, no unknown (N=2507)")

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
data.model = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry %in% c("AFR", "AMR", "EUR", "EAS", "SAS")),]

ASD_effect[[7]]=PGS.TRI(pgs_offspring = data.model$pgs_child_PC_sd,
                        pgs_mother = data.model$pgs_mother_PC_sd,
                        pgs_father = data.model$pgs_father_PC_sd, 
                        GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
res_direct_PC[7,] = ASD_effect[[7]]$res_beta
res_indirect_PC[7,] = ASD_effect[[7]]$res_delta

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


direct_effect_heter = rbind(res_direct,res_direct_PC,res_direct_PC_mean,res_direct_1KG)
indirect_effect_heter = rbind(res_indirect,res_indirect_PC,res_indirect_PC_mean,res_indirect_1KG)


save(direct_effect_all,direct_effect_homo,indirect_effect_all,indirect_effect_homo,direct_effect_heter,indirect_effect_heter,file="C:/Users/ziqia/Desktop/work/family/data/autism/v3/results_ASD_all.RData")
write.csv(direct_effect_all,file="C:/Users/ziqia/Desktop/work/family/data/autism/v3/results_ASD_direct_all.csv")
write.csv(direct_effect_homo,file="C:/Users/ziqia/Desktop/work/family/data/autism/v3/results_ASD_direct_homo.csv")
write.csv(direct_effect_heter,file="C:/Users/ziqia/Desktop/work/family/data/autism/v3/results_ASD_direct_heter.csv")
write.csv(indirect_effect_all,file="C:/Users/ziqia/Desktop/work/family/data/autism/v3/results_ASD_indirect_all.csv")
write.csv(indirect_effect_homo,file="C:/Users/ziqia/Desktop/work/family/data/autism/v3/results_ASD_indirect_homo.csv")
write.csv(indirect_effect_heter,file="C:/Users/ziqia/Desktop/work/family/data/autism/v3/results_ASD_indirect_heter.csv")

