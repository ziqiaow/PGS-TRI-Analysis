library(data.table)

parental_geno = fread("/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/data/chr1.raw")
for(i in 2:22){
  
  tmp = fread(paste0("/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/data/chr",i,".raw"))
  identical(tmp$IID,parental_geno$IID)
  parental_geno = cbind(parental_geno,tmp[,-c(1:6)])
  

}
save(parental_geno,file="/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/data/matched_parents.RData")
fwrite(parental_geno,file="/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/data/matched_parents.raw",row.names=TRUE, quote=F,col.names=T)



#Next clean up weight files
sum_stats = fread(paste0("/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/sum_stats_chr",1,".txt"))
snps = fread(paste0("/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/clumped_snps_chr",1,".txt"),sep = "\t")
sum_stats = sum_stats[match(snps$x,sum_stats$ID),]
sum_stats$updated_weight = sum_stats$BETA * (-1)
for(i in 2:22){
  
  sum_stats_tmp = fread(paste0("/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/sum_stats_chr",i,".txt"))
  snps = fread(paste0("/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/clumped_snps_chr",i,".txt"),sep = "\t")
 
  sum_stats_tmp = sum_stats_tmp[match(snps$x,sum_stats_tmp$ID),]
  sum_stats_tmp$updated_weight = sum_stats_tmp$BETA * (-1)
  
  sum_stats = rbind(sum_stats,sum_stats_tmp)
  

  
}
#match the order with parental geno
id = colnames(parental_geno)[-c(1:6)]
sum_stats$ID_ref = paste0(sum_stats$ID,"_",sum_stats$REF)
sum_stats = sum_stats[match(id,sum_stats$ID_ref),]
#save this file for calculating scores in plink
write.table(sum_stats,file=paste0("/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/weight_final.txt"),row.names = F,col.names = T,quote = F)
save(sum_stats,file="/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/weight_final.RData")


#Check the parental genotype and children's allele frequency
parental_geno = load("/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/data/matched_parents.RData")
#parental_geno = parental_geno[,-1]
#too slow to re-read into the environment, directly use it
load("/dcs04/nilanjan/data/zwang/family/simulation/matched_parents_pheno.RData")
id_parental$fam = 1:145956

#colnames(parental_geno)[-c(1:6)] = paste0("snp",1:4483)
father_geno = data.frame(parental_geno[match(id_parental$X1,parental_geno$IID),-c(1:6)])
mother_geno = data.frame(parental_geno[match(id_parental$X2,parental_geno$IID),-c(1:6)])
#rm(parental_geno)
save(father_geno,mother_geno,id_parental,file="/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/parental_geno_pfile.RData")


#load weights
weight = fread("/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/weight_final.txt")

#calculate the parental PGS using R instead so that we can match the algorithms between parents and children
#father_geno = father_geno[,match(weight$snpindex,colnames(father_geno))]
#mother_geno = mother_geno[,match(weight$snpindex,colnames(mother_geno))]


#calculate the parental PGS using R instead so that we can match the algorithms between parents and children
father_geno=apply(father_geno,2, function(x) { 
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

mother_geno=apply(mother_geno,2, function(x) { 
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

pgs_f = as.numeric(data.matrix(father_geno) %*% weight$updated_weight)
pgs_m = as.numeric(data.matrix(mother_geno) %*% weight$updated_weight)

names(pgs_f) = id_parental$X1
names(pgs_m) = id_parental$X2
save(pgs_f,pgs_m,id_parental,file="/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/pgs/PGS_parents_id_R_pfile.RData")




#for match_geo parental data

#Check the parental genotype and children's allele frequency
#parental_geno = fread("/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440/plink/data_pfile/matched_parents.raw")
#too slow to re-read into the environment, directly use it
load("/dcs04/nilanjan/data/zwang/family/simulation/match_geo/matched_parents_pheno.RData")
id_parental$fam = 1:150253

#colnames(parental_geno)[-c(1:6)] = paste0("snp",1:19452)
father_geno = data.frame(parental_geno[match(id_parental$X1,parental_geno$IID),-c(1:6)])
mother_geno = data.frame(parental_geno[match(id_parental$X2,parental_geno$IID),-c(1:6)])
#rm(parental_geno)
save(father_geno,mother_geno,id_parental,file="/dcs04/nilanjan/data/zwang/family/simulation/match_geo/new_score/clump/pgs/parental_geno_pfile.RData")


#load weights
weight = fread("/dcs04/nilanjan/data/zwang/family/simulation/new_score/PGS002440_clump/weight_final.txt")

#calculate the parental PGS using R instead so that we can match the algorithms between parents and children
# father_geno = father_geno[,match(weight$snpindex,colnames(father_geno))]
# mother_geno = mother_geno[,match(weight$snpindex,colnames(mother_geno))]



#calculate the parental PGS using R instead so that we can match the algorithms between parents and children
father_geno=apply(father_geno,2, function(x) { 
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

mother_geno=apply(mother_geno,2, function(x) { 
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
})

pgs_f = as.numeric(data.matrix(father_geno) %*% weight$updated_weight)
pgs_m = as.numeric(data.matrix(mother_geno) %*% weight$updated_weight)

names(pgs_f) = id_parental$X1
names(pgs_m) = id_parental$X2
save(pgs_f,pgs_m,id_parental,file="/dcs04/nilanjan/data/zwang/family/simulation/match_geo/new_score/clump/pgs/PGS_parents_id_R_pfile.RData")





