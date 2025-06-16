#First, simulate 100,000 families based on 20 generations of assortative mating
#Program: https://snipar.readthedocs.io/en/latest/simulation.html
#Reference: https://www.nature.com/articles/s41588-022-01085-0
#location on JHPCE:
#/dcs04/nilanjan/data/zwang/family/simulation/snipar/simulated_data
#code:
#/dcs04/nilanjan/data/zwang/family/simulation/snipar/simulated_data/simulated_data.sh
#simulate.py 1000 0.5 /dcs04/nilanjan/data/zwang/family/simulation/snipar/simulated_data/ --nfam 100000 --n_am 20 --r_par 0.5 --save_par_gts

#next, use causal_effects.txt to calculate PGS in simulated genotypes of parents and children

#finally, export these PGS and simulate disease status
library(data.table)
setwd("~/work/family/simulation/snipar")
pedigree = fread("pedigree.txt")
pgs_children = fread("PGS_children.sscore")
pgs_parents = fread("PGS_parents.sscore")
even_indexes<-seq(2,200000,2)
odd_indexes<-seq(1,199999,2)

pgs_children_kid1 = pgs_children[odd_indexes,]
pgs_children_kid2 = pgs_children[even_indexes,]
hist(pgs_children_kid1$SCORE1_SUM)
hist(pgs_children_kid2$SCORE1_SUM)
summary(pgs_children_kid1$SCORE1_SUM)
summary(pgs_children_kid2$SCORE1_SUM)


pgs_m_id = pedigree$MOTHER_ID[match(pgs_children_kid1$IID,pedigree$IID)]
pgs_m = pgs_parents[match(pgs_m_id,pgs_parents$IID),]
pgs_f_id = pedigree$FATHER_ID[match(pgs_children_kid1$IID,pedigree$IID)]
pgs_f = pgs_parents[match(pgs_f_id,pgs_parents$IID),]


hist(pgs_m$SCORE1_SUM)
hist(pgs_f$SCORE1_SUM)
summary(pgs_m$SCORE1_SUM)
summary(pgs_f$SCORE1_SUM)

#identical(pgs_m$`#FID`, pgs_f$`#FID`)
cor(pgs_m$SCORE1_SUM,pgs_f$SCORE1_SUM)
cor(pgs_children_kid1$SCORE1_SUM,pgs_f$SCORE1_SUM)
cor(pgs_children_kid1$SCORE1_SUM,pgs_m$SCORE1_SUM)


#export parental phenotype
pheno_parent = 0.5*(pedigree$PHENO[match(pgs_m$IID,pedigree$IID)] + pedigree$PHENO[match(pgs_f$IID,pedigree$IID)])



save(pgs_m,pgs_f,pgs_children_kid1,pheno_parent,file="PGS_family.RData")



#finally, export these PGS and simulate disease status
library(data.table)
setwd("~/work/family/simulation/snipar")
pgs_children = fread("PGS_children_v1.sscore")
pgs_parents = fread("PGS_parents_v1.sscore")
even_indexes<-seq(2,200000,2)
odd_indexes<-seq(1,199999,2)

pgs_children_kid1 = pgs_children[odd_indexes,]
pgs_children_kid2 = pgs_children[even_indexes,]
hist(pgs_children_kid1$SCORE1_SUM)
hist(pgs_children_kid2$SCORE1_SUM)
summary(pgs_children_kid1$SCORE1_SUM)
summary(pgs_children_kid2$SCORE1_SUM)


pgs_m_id = pedigree$MOTHER_ID[match(pgs_children_kid1$IID,pedigree$IID)]
pgs_m = pgs_parents[match(pgs_m_id,pgs_parents$IID),]
pgs_f_id = pedigree$FATHER_ID[match(pgs_children_kid1$IID,pedigree$IID)]
pgs_f = pgs_parents[match(pgs_f_id,pgs_parents$IID),]


hist(pgs_m$SCORE1_SUM)
hist(pgs_f$SCORE1_SUM)
summary(pgs_m$SCORE1_SUM)
summary(pgs_f$SCORE1_SUM)

#identical(pgs_m$`#FID`, pgs_f$`#FID`)
cor(pgs_m$SCORE1_SUM,pgs_f$SCORE1_SUM)
cor(pgs_children_kid1$SCORE1_SUM,pgs_f$SCORE1_SUM)
cor(pgs_children_kid1$SCORE1_SUM,pgs_m$SCORE1_SUM)

save(pgs_m,pgs_f,pgs_children_kid1,file="PGS_family_v1.RData")
