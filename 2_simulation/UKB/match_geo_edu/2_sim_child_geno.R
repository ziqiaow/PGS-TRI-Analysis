#simulate children's genotype based on parental genotype
#May 5, 2024
#------------------------------------------------------------
#Check the parental genotype and children's allele frequency
path <- "/dcs04/nilanjan/data/zwang/family/simulation/ukb_05052024/child_geno/"

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
task_id <- as.numeric(slurm_arrayid)

print(task_id)
set.seed(task_id)

library(data.table)
parental_geno = fread("/dcs04/nilanjan/data/zwang/family/simulation/pgs/data/matched_parents_pruned.raw")
load("/dcs04/nilanjan/data/zwang/family/simulation/matched_parents_pheno.RData")
id_parental$fam = 1:145956

father_geno = data.frame(parental_geno[match(id_parental$X1,parental_geno$IID),-c(1,2,3,4,5,6)])
mother_geno = data.frame(parental_geno[match(id_parental$X2,parental_geno$IID),-c(1,2,3,4,5,6)])
rm(parental_geno)

myfunc = function(x,y){
  if((is.na(x) | is.na(y))){return(NA)}
  if(x==y & x!=1){
    
    return(x)
    
  } else if(abs(x-y)==2){
    
    return(1) 
    
  }else if(abs(x-y)==1){
    
    return( (rbinom(1,1,0.5)+min(x,y)) )
    
  } else {
    
    return(sample(c(0, 1, 2), size = 1, prob = c(0.25, 0.5, 0.25)))
  }
}


myfunc2 = function(a,b) {mapply(myfunc, a, b)}
child_geno = mapply(function(x,y) myfunc2(x,y), as.data.frame(father_geno), as.data.frame(mother_geno))
#about 10 min 

save(child_geno,file=paste0(path,"child_geno_",task_id,".RData"))




weight = fread("/dcs04/nilanjan/data/zwang/family/simulation/pgs/ea_score_pruned.txt")
variants = fread("/dcs04/nilanjan/data/zwang/family/simulation/pgs/ea_score_parents_pruned.sscore.vars",header = F)
weight = weight[match(variants$V1,weight$V1),]

#first impute missing genotype data for children
child_geno=apply(child_geno,2, function(x) { 
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

pgs_c = as.numeric(child_geno %*% weight$V3)
save(pgs_c,file=paste0(path,"child_pgs_",task_id,".RData"))


