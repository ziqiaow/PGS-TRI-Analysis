#simulation work for case-parent trio
#Prospectively simulate the data
#Aug 12, 2025
#---------------------------------------------

path<-"/dcs04/nilanjan/data/zwang/family/simulation/match_geo/new_score/clump/PS/PGS_main/"
source("/users/zwang4/family/R/PGS-TRI.R")
source("/users/zwang4/family/R/pTDT.R")
source("/users/zwang4/family/R/sim_disease_risk_UKB.R")

load("/dcs04/nilanjan/data/zwang/family/simulation/match_geo/matched_parents_pheno.RData")
load("/dcs04/nilanjan/data/zwang/family/simulation/match_geo/new_score/clump/pgs/PGS_parents_id_R_pfile.RData")

sigma = 0.5*(pgs_m - pgs_f)^2

# Calculate quartiles and IQR
Q1 <- quantile(sigma, 0.25)
Q3 <- quantile(sigma, 0.75)
IQR_val <- Q3 - Q1

# Define outlier bounds
lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

# Remove outliers
id <- which(sigma >= lower_bound & sigma <= upper_bound)
pheno_father = pheno_father[id,]
pheno_mother = pheno_mother[id,]
id_parental = id_parental[id,]
pgs_f = pgs_f[id]
pgs_m = pgs_m[id]


bmi = pheno_mother$BMI
means_bmi = aggregate(pheno_mother$BMI,by = list(pheno_mother$education_match),mean,na.rm=T)
bmi[which(is.na(bmi))] = means_bmi$x[match(as.character(pheno_mother$education_match[which(is.na(bmi))]),as.character(means_bmi$Group.1))]
bmi=scale(bmi)
bmi = bmi*(-1)


pheno_child = pheno_father[,c(2,3,4)]
pheno_child$birth.x = scale(pheno_child$birth.x)
pheno_child$birth.y = scale(pheno_child$birth.y)
pheno_child$assessment_center = factor(pheno_child$assessment_center)

pc_child = (pheno_mother[,11:20] + pheno_father[,11:20])/2

#This part does not change
sim_N=1000
fam_size1=1000
res_fam_size1=res_glm_fam_mf=res_glm_fam_size1=res_ptdt_size1=res_glm_pc10_fam_size1=res_glm_pc10_geo_fam_size1=list()
res_pc_father=res_pc_mother=res_pc_famo=list()


h=0
for(i in 1:125){
 #set.seed(i)
  #load simulated children's pgs
  load(paste0("/dcs04/nilanjan/data/zwang/family/simulation/match_geo/new_score/clump/child_geno/child_pgs_",i,".RData"))
  pgs_c = pgs_c[id]
  
  for(j in 1:8){
  h=h+1
  
  #first without PCs in the intercept term
  dat = sim_prospective_population(pgs_c = pgs_c , pgs_f = pgs_f,pgs_m = pgs_m,e_strat = bmi,cor_strat = 2,alpha_fam = -6.3,pc=F,mu_pc = mu_sim,betaG_normPRS=0)
  
  #First calculate family size 1 = 1000
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size1,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size1,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(dat$pgs_fam[id1,1])),as.numeric(as.vector(dat$pgs_fam[id2,1])))
  PRS_fam=dat$pgs_fam[id2,]
  pheno_sim = pheno_child[c(id1,id2),]
  
  G_mother = c(as.numeric(as.vector(dat$pgs_fam[id1,2])),as.numeric(as.vector(dat$pgs_fam[id2,2])))
  G_father = c(as.numeric(as.vector(dat$pgs_fam[id1,3])),as.numeric(as.vector(dat$pgs_fam[id2,3])))
  
  res_fam_size1[[h]]=PGS.TRI(pgs_offspring = PRS_fam[,1], pgs_mother = PRS_fam[,2], pgs_father=PRS_fam[,3], side=2,GxE_int = FALSE)$Coefficients_direct
  
  fit <- glm(D ~ G_logistic+G_mother+G_father, family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_fam_mf[[h]] <- summary(fit)$coef
  
   fit <- glm(D ~ G_logistic, family = binomial(), model = FALSE, x = FALSE, y = FALSE)
   res_glm_fam_size1[[h]]=summary(fit)$coef
  
  pc_sim = pc_child[c(id1,id2),1:10]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc10 <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_pc10_fam_size1[[h]]=summary(fit_pc10)$coef
  

  dat_logit = data.frame(cbind(G_logistic,pc_sim,pheno_sim))
  fit_pc10_geo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_pc10_geo_fam_size1[[h]]=summary(fit_pc10_geo)$coef
  
  
  pc_sim = pheno_mother[c(id1,id2),11:20]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_mo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_mother[[h]]=summary(fit_pc_mo)$coef
  
  pc_sim = pheno_father[c(id1,id2),11:20]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_fa <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_father[[h]]=summary(fit_pc_fa)$coef
  
  pc_sim = cbind(pheno_father[c(id1,id2),11:20],pheno_mother[c(id1,id2),11:20])
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_famo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_famo[[h]]=summary(fit_pc_famo)$coef
  
  # dat_logit = data.frame(cbind(G_logistic,edu_sim))
  # fit_edu <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  # res_glm_edu_fam_size1[[h]]=summary(fit_edu)$coef

  
  res_ptdt_size1[[h]]=ptdt(PRS_fam[,1],PRS_fam[,2], PRS_fam[,3],side0=2)$statistic
  
  
  }
}

save(res_fam_size1,res_glm_fam_size1,res_ptdt_size1,res_pc_father,res_pc_mother,res_glm_fam_mf,
     res_glm_pc10_fam_size1,res_glm_pc10_geo_fam_size1,res_pc_famo,file=paste0(path,"res_beta0_size1k_trim_multichild",".rda"))





#This part does not change
sim_N=1000
fam_size1=2000
res_fam_size1=res_glm_fam_mf=res_glm_fam_size1=res_ptdt_size1=res_glm_pc10_fam_size1=res_glm_pc10_geo_fam_size1=list()
res_pc_father=res_pc_mother=res_pc_famo=list()
#increase the sample size
h=0
for(i in 1:125){
  #set.seed(i)
  #load simulated children's pgs
  load(paste0("/dcs04/nilanjan/data/zwang/family/simulation/match_geo/new_score/clump/child_geno/child_pgs_",i,".RData"))
  pgs_c = pgs_c[id]
  
  for(j in 1:8){
    h=h+1

  #first without PCs in the intercept term
  dat = sim_prospective_population(pgs_c = pgs_c , pgs_f = pgs_f,pgs_m = pgs_m,e_strat = bmi,cor_strat = 2,alpha_fam = -5.6,pc=F,mu_pc = mu_sim,betaG_normPRS=0)
  
  #First calculate family size 1 = 1000
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size1,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size1,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(dat$pgs_fam[id1,1])),as.numeric(as.vector(dat$pgs_fam[id2,1])))
  #G_logistic = scale(G_logistic)
  PRS_fam=dat$pgs_fam[id2,]
  # tmp_all = c(PRS_fam[,1],PRS_fam[,2],PRS_fam[,3])
  # PRS_fam = apply(PRS_fam, 2, function(x) (x - mean(tmp_all))/sd(tmp_all))
  # edu_sim = edu[c(id1,id2)]
  pheno_sim = pheno_child[c(id1,id2),]
  G_mother = c(as.numeric(as.vector(dat$pgs_fam[id1,2])),as.numeric(as.vector(dat$pgs_fam[id2,2])))
  G_father = c(as.numeric(as.vector(dat$pgs_fam[id1,3])),as.numeric(as.vector(dat$pgs_fam[id2,3])))
  
  res_fam_size1[[h]]=PGS.TRI(pgs_offspring = PRS_fam[,1], pgs_mother = PRS_fam[,2], pgs_father=PRS_fam[,3], side=2,GxE_int = FALSE)$Coefficients_direct
  
  fit <- glm(D ~ G_logistic+G_mother+G_father, family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_fam_mf[[h]] <- summary(fit)$coef
  
  fit <- glm(D ~ G_logistic, family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_fam_size1[[h]]=summary(fit)$coef
  
  pc_sim = pc_child[c(id1,id2),1:10]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc10 <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_pc10_fam_size1[[h]]=summary(fit_pc10)$coef
  
  
  dat_logit = data.frame(cbind(G_logistic,pc_sim,pheno_sim))
  fit_pc10_geo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_pc10_geo_fam_size1[[h]]=summary(fit_pc10_geo)$coef
  
  
  pc_sim = pheno_mother[c(id1,id2),11:20]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_mo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_mother[[h]]=summary(fit_pc_mo)$coef
  
  pc_sim = pheno_father[c(id1,id2),11:20]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_fa <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_father[[h]]=summary(fit_pc_fa)$coef
  
  pc_sim = cbind(pheno_father[c(id1,id2),11:20],pheno_mother[c(id1,id2),11:20])
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_famo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_famo[[h]]=summary(fit_pc_famo)$coef
  
  # dat_logit = data.frame(cbind(G_logistic,edu_sim))
  # fit_edu <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  # res_glm_edu_fam_size1[[h]]=summary(fit_edu)$coef
  
  
  res_ptdt_size1[[h]]=ptdt(PRS_fam[,1],PRS_fam[,2], PRS_fam[,3],side0=2)$statistic
  
  
 }
}

save(res_fam_size1,res_glm_fam_size1,res_ptdt_size1,res_pc_father,res_pc_mother,res_glm_fam_mf,
     res_glm_pc10_fam_size1,res_glm_pc10_geo_fam_size1,res_pc_famo,file=paste0(path,"res_beta0_size2k_trim_multichild",".rda"))



