#Check parental PGS correlation and how well it predicts education
#UKB
#Aug 5, 2025
#---------------------------------------------
#supplemental table 

base::load("~/work/family/simulation/UKB/new_score/clump/pgs/PGS_parents_id_R_pfile.RData")
base::load("~/work/family/simulation/UKB/matched_parents_pheno.RData")

pheno_all = rbind(pheno_mother,pheno_father)
pgs_all = c(pgs_m,pgs_f)
cor.test(pheno_all$education,pgs_all) #-0.12302  
test=lm(pheno_all$education ~ pgs_all)
summary(test) #R-squared:  0.01513    


cor.test(pheno_mother$education,pgs_m) #-0.130325
test=lm(pheno_mother$education ~ pgs_m)
summary(test) #R-squared:  0.01698  

cor.test(pheno_father$education,pgs_f) #-0.1154453
test=lm(pheno_father$education ~ pgs_f)
summary(test) #R-squared:  0.01332  
cor.test(pgs_f,pgs_m) #0.1357339  

load("~/work/family/simulation/UKB/new_score/clump/child_geno/child_pgs_1.RData")
cor.test(pgs_c,pgs_m)
cor.test(pgs_c,pgs_f)

