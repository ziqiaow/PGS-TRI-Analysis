#simulate children's genotype based on parental genotype
#Aug 10, 2025
#------------------------------------------------------------
#Check the parental genotype and children's allele frequency
path <- "/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/child_geno/"

library(data.table)
load("/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/pgs/parental_geno_pfile.RData")

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
task_id <- as.numeric(slurm_arrayid)

print(task_id)
set.seed(task_id)


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




#calculate child PGS
#load weights
weight = fread("/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/weight_final.txt")

child_geno=apply(child_geno,2, function(x) { 
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

pgs_c = as.numeric(data.matrix(child_geno) %*% weight$updated_weight)
save(pgs_c,file=paste0(path,"child_pgs_",task_id,".RData"))
save(child_geno,file=paste0(path,"child_geno_",task_id,".RData"))



