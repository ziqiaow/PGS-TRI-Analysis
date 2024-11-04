library(data.table)
#Check the parental genotype and children's allele frequency
parental_geno = fread("/dcs04/nilanjan/data/zwang/family/simulation/match_geo/plink/data/matched_parents_pruned.raw")
load("/dcs04/nilanjan/data/zwang/family/simulation/match_geo/matched_parents_pheno.RData")
id_parental$fam = 1:150253

father_geno = data.frame(parental_geno[match(id_parental$X1,parental_geno$IID),-c(1,2,3,4,5,6)])
mother_geno = data.frame(parental_geno[match(id_parental$X2,parental_geno$IID),-c(1,2,3,4,5,6)])
rm(parental_geno)


#load weights
weight = fread("/dcs04/nilanjan/data/zwang/family/simulation/pgs/ea_score_pruned.txt")
variants = fread("/dcs04/nilanjan/data/zwang/family/simulation/pgs/ea_score_parents_pruned.sscore.vars",header = F)
weight = weight[match(variants$V1,weight$V1),]

#calculate the parental PGS using R instead so that we can match the algorithms between parents and children
father_geno=apply(father_geno,2, function(x) { 
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})
pgs_f = as.numeric(father_geno %*% weight$V3)

mother_geno=apply(mother_geno,2, function(x) { 
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})
pgs_m = as.numeric(mother_geno %*% weight$V3)

names(pgs_f) = id_parental$X1
names(pgs_m) = id_parental$X2
save(pgs_f,pgs_m,id_parental,file="/dcs04/nilanjan/data/zwang/family/simulation/match_geo/PGS_parents_id_R.RData")

