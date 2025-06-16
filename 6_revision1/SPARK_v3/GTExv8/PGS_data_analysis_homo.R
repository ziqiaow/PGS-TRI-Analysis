#PGS analysis for omics
library(data.table)
gtex_list = fread("/dcs04/nilanjan/data/zwang/family/gtex/data/Tissue_list_GTex_V8_analyze.txt",header = F)
data_name = gtex_list$V1


# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
task_id <- as.numeric(slurm_arrayid)


source("/users/zwang4/family/R/PGS-TRI.R")
source("/users/zwang4/family/R/pTDT.R")

#load family information
load("/dcs04/nilanjan/data/zwang/family/autism/roles_id_complete_casetrio_ancestry.RData")
id_final = id_final[,c(1,9,10,38:41)]
id_final = id_final[which(id_final$ancestry == id_final$ancestry_f & id_final$ancestry == id_final$ancestry_m),]
dim(id_final)

#load omics PGS data
load(paste0("/dcs04/nilanjan/data/zwang/family/gtex/data/results/",data_name[task_id],".RData"))
p = dim(omics_pgs)[2] - 4
direct_effect=array(0,c(p,3*6))
indirect_effect=array(0,c(p,3*6))
colnames(direct_effect)=colnames(indirect_effect)=c("beta_all","se_all","p_all",
                                                    "beta_EUR","se_EUR","p_EUR",
                                                    "beta_AFR","se_AFR","p_AFR",
                                                    "beta_AMR","se_AMR","p_AMR",
                                                    "beta_EAS","se_EAS","p_EAS",
                                                    "beta_SAS","se_SAS","p_SAS")

group = c("EUR","AFR","AMR","EAS","SAS")

for(i in 1:p){
  print(i)
  
  id_final$pgs_child = omics_pgs[match(id_final$subject_sp_id,omics_pgs$IID),i+4]
  id_final$pgs_mother = omics_pgs[match(id_final$biomother_sp_id,omics_pgs$IID),i+4]
  id_final$pgs_father = omics_pgs[match(id_final$biofather_sp_id,omics_pgs$IID),i+4]

  #Analysis for all population (including unknown), no standardization
  prs_pgscatalog_complete = id_final
  
  #all
  res = PGS.TRI(pgs_offspring = prs_pgscatalog_complete$pgs_child,
                   pgs_mother = prs_pgscatalog_complete$pgs_mother,
                   pgs_father = prs_pgscatalog_complete$pgs_father, 
                   GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
  direct_effect[i,1:3]=res$res_beta[c(1,2,4)]
  indirect_effect[i,1:3]=res$res_delta[c(1,2,4)]
  
  
  for(j in 2:6){
    
    data.model = prs_pgscatalog_complete[which(prs_pgscatalog_complete$ancestry == group[j-1]),]
    res = PGS.TRI(pgs_offspring = data.model$pgs_child,
                            pgs_mother = data.model$pgs_mother,
                            pgs_father = data.model$pgs_father, 
                            GxE_int = FALSE, smalltriosize = F, parental_indirect = T)
    direct_effect[i,(j*3-2):(j*3)] = res$res_beta[c(1,2,4)]
    indirect_effect[i,(j*3-2):(j*3)] = res$res_delta[c(1,2,4)]
    
    
  }
}

rownames(direct_effect)=rownames(indirect_effect)=colnames(omics_pgs)[-c(1:4)] 
save(direct_effect,indirect_effect,file=paste0("/dcs04/nilanjan/data/zwang/family/gtex/data/results/",data_name[task_id],"_results_homo.RData"))
write.csv(direct_effect,file=paste0("/dcs04/nilanjan/data/zwang/family/gtex/data/results/",data_name[task_id],"_direct_homo.csv"))
write.csv(indirect_effect,file=paste0("/dcs04/nilanjan/data/zwang/family/gtex/data/results/",data_name[task_id],"_indirect_homo.csv"))