#Summarize the PGS scores for each omic type
library(data.table)

gtex_list = fread("/dcs04/nilanjan/data/zwang/family/gtex/data/Tissue_list_GTex_V8_analyze.txt",header = F)
data_name = gtex_list$V1

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
task_id <- as.numeric(slurm_arrayid)


#load("/dcs04/nilanjan/data/zwang/family/autism/roles_id_complete_casetrio.RData")
# dat <- readRDS("/dcs04/nilanjan/data/zwang/family/autism/plink/pgs/redo/pgs_ASD_final_ancestry.rds")
# 
# id_final$ancestry = dat$superclass[match(id_final$subject_sp_id,dat$IID)]
# id_final$ancestry_m = dat$superclass[match(id_final$biomother_sp_id,dat$IID)]
# id_final$ancestry_f = dat$superclass[match(id_final$biofather_sp_id,dat$IID)]
# id_final$ancestry_family = paste0(id_final$ancestry,",m=",id_final$ancestry_m,",f=",id_final$ancestry_f)
# save(id_final,file="/dcs04/nilanjan/data/zwang/family/autism/roles_id_complete_casetrio_ancestry.RData")
# 

#for(type in gtex_list$V1){

  #load trio information
  load("/dcs04/nilanjan/data/zwang/family/autism/roles_id_complete_casetrio_ancestry.RData")
  
id_final = id_final[,c(1,9,10,38:41)]
id = c(id_final$subject_sp_id,id_final$biomother_sp_id,id_final$biofather_sp_id)
filenames <- list.files(paste0("/dcs04/nilanjan/data/zwang/family/gtex/data/scores/",data_name[task_id]), pattern=".sscore$", full.names=TRUE)


omics_pgs = data.frame(fread(filenames[1]))
omics_pgs = omics_pgs[match(id,omics_pgs$IID),]

for(i in 2:length(filenames)){
  tmp = data.frame(fread(filenames[i]))
  tmp = tmp[match(id,tmp$IID),]
  omics_pgs = cbind(omics_pgs,tmp$SCORE1_SUM)
}

name_tmp=gsub(".*/","",filenames)
name_tmp=gsub("\\..*","",name_tmp)
colnames(omics_pgs)[-c(1:4)] = name_tmp
save(omics_pgs,file=paste0("/dcs04/nilanjan/data/zwang/family/gtex/data/results/",data_name[task_id],".RData"))
#}


