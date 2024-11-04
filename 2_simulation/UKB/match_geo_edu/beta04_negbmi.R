#simulation work for case-parent trio
#Prospectively simulate the data
#May 5, 2024
#---------------------------------------------

path<-"/dcs04/nilanjan/data/zwang/family/simulation/ukb_05132024/"
source("/users/zwang4/family/R/PGS-TRI.R")
source("/users/zwang4/family/R/sim_disease_risk_UKB.R")
load("/dcs04/nilanjan/data/zwang/family/simulation/matched_parents_pheno.RData")
load("/dcs04/nilanjan/data/zwang/family/simulation/PGS_parents_id_R.RData")

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
res_fam_size1=res_glm_fam_size1=res_ptdt_size1=res_glm_pc10_fam_size1=res_glm_pc10_geo_fam_size1=list()
res_pc_father=res_pc_mother=res_pc_famo=list()





for(i in 1:sim_N){
 set.seed(i)
  #load simulated children's pgs
  load(paste0("/dcs04/nilanjan/data/zwang/family/simulation/ukb_05052024/child_geno/child_pgs_",i,".RData"))

  
  #first without PCs in the intercept term
  dat = sim_prospective_population(pgs_c = pgs_c , pgs_f = pgs_f,pgs_m = pgs_m,e_strat = bmi,cor_strat = 2,alpha_fam = -6.3,pc=F,mu_pc = mu_sim,betaG_normPRS=0.4)
  
  #First calculate family size 1 = 1000
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size1,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size1,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(dat$pgs_fam[id1,1])),as.numeric(as.vector(dat$pgs_fam[id2,1])))
  PRS_fam=dat$pgs_fam[id2,]
  pheno_sim = pheno_child[c(id1,id2),]
  
   res_fam_size1[[i]]=PGS.TRI(pgs_offspring = PRS_fam[,1], pgs_mother = PRS_fam[,2], pgs_father=PRS_fam[,3], side=2,GxE_int = FALSE)$res_beta
   fit <- glm(D ~ G_logistic, family = binomial(), model = FALSE, x = FALSE, y = FALSE)
   res_glm_fam_size1[[i]]=summary(fit)$coef
  
  pc_sim = pc_child[c(id1,id2),1:10]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc10 <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_pc10_fam_size1[[i]]=summary(fit_pc10)$coef
  

  dat_logit = data.frame(cbind(G_logistic,pc_sim,pheno_sim))
  fit_pc10_geo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_pc10_geo_fam_size1[[i]]=summary(fit_pc10_geo)$coef
  
  
  pc_sim = pheno_mother[c(id1,id2),11:20]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_mo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_mother[[i]]=summary(fit_pc_mo)$coef
  
  pc_sim = pheno_father[c(id1,id2),11:20]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_fa <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_father[[i]]=summary(fit_pc_fa)$coef
  
  pc_sim = cbind(pheno_father[c(id1,id2),11:20],pheno_mother[c(id1,id2),11:20])
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_famo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_famo[[i]]=summary(fit_pc_famo)$coef
  

  res_ptdt_size1[[i]]=ptdt(PRS_fam[,1],PRS_fam[,2], PRS_fam[,3],side0=2)$res_beta
  
  
}

save(res_fam_size1,res_glm_fam_size1,res_ptdt_size1,res_pc_father,res_pc_mother,
     res_glm_pc10_fam_size1,res_glm_pc10_geo_fam_size1,res_pc_famo,file=paste0(path,"res_beta04_size1k",".rda"))





#This part does not change
sim_N=1000
fam_size1=2000
res_fam_size1=res_glm_fam_size1=res_ptdt_size1=res_glm_pc10_fam_size1=res_glm_pc10_geo_fam_size1=list()
res_pc_father=res_pc_mother=res_pc_famo=list()
#increase the sample size
for(i in 1:sim_N){
  set.seed(i)
  #load simulated children's pgs
  load(paste0("/dcs04/nilanjan/data/zwang/family/simulation/ukb_05052024/child_geno/child_pgs_",i,".RData"))

  #first without PCs in the intercept term
  dat = sim_prospective_population(pgs_c = pgs_c , pgs_f = pgs_f,pgs_m = pgs_m,e_strat = bmi,cor_strat = 2,alpha_fam = -5.5,pc=F,mu_pc = mu_sim,betaG_normPRS=0.4)
  
  #First calculate family size 1 = 1000
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size1,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size1,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(dat$pgs_fam[id1,1])),as.numeric(as.vector(dat$pgs_fam[id2,1])))
  PRS_fam=dat$pgs_fam[id2,]
  pheno_sim = pheno_child[c(id1,id2),]
  
  res_fam_size1[[i]]=PGS.TRI(pgs_offspring = PRS_fam[,1], pgs_mother = PRS_fam[,2], pgs_father=PRS_fam[,3], side=2,GxE_int = FALSE)$res_beta
  fit <- glm(D ~ G_logistic, family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_fam_size1[[i]]=summary(fit)$coef
  
  pc_sim = pc_child[c(id1,id2),1:10]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc10 <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_pc10_fam_size1[[i]]=summary(fit_pc10)$coef
  
  
  dat_logit = data.frame(cbind(G_logistic,pc_sim,pheno_sim))
  fit_pc10_geo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_pc10_geo_fam_size1[[i]]=summary(fit_pc10_geo)$coef
  
  
  pc_sim = pheno_mother[c(id1,id2),11:20]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_mo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_mother[[i]]=summary(fit_pc_mo)$coef
  
  pc_sim = pheno_father[c(id1,id2),11:20]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_fa <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_father[[i]]=summary(fit_pc_fa)$coef
  
  pc_sim = cbind(pheno_father[c(id1,id2),11:20],pheno_mother[c(id1,id2),11:20])
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_famo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_famo[[i]]=summary(fit_pc_famo)$coef
  
  
  res_ptdt_size1[[i]]=ptdt(PRS_fam[,1],PRS_fam[,2], PRS_fam[,3],side0=2)$res_beta
  
  
}

save(res_fam_size1,res_glm_fam_size1,res_ptdt_size1,res_pc_father,res_pc_mother,
     res_glm_pc10_fam_size1,res_glm_pc10_geo_fam_size1,res_pc_famo,file=paste0(path,"res_beta04_size2k",".rda"))



#This part does not change
sim_N=1000
fam_size1=3000
res_fam_size1=res_glm_fam_size1=res_ptdt_size1=res_glm_pc10_fam_size1=res_glm_pc10_geo_fam_size1=list()
res_pc_father=res_pc_mother=res_pc_famo=list()
#increase the sample size
for(i in 1:sim_N){
  set.seed(i)
  #load simulated children's pgs
  load(paste0("/dcs04/nilanjan/data/zwang/family/simulation/ukb_05052024/child_geno/child_pgs_",i,".RData"))

  #first without PCs in the intercept term
  dat = sim_prospective_population(pgs_c = pgs_c , pgs_f = pgs_f,pgs_m = pgs_m,e_strat = bmi,cor_strat = 2,alpha_fam = -5.2,pc=F,mu_pc = mu_sim,betaG_normPRS=0.4)
  
  #First calculate family size 1 = 1000
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size1,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size1,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(dat$pgs_fam[id1,1])),as.numeric(as.vector(dat$pgs_fam[id2,1])))
  PRS_fam=dat$pgs_fam[id2,]
  pheno_sim = pheno_child[c(id1,id2),]
  
  res_fam_size1[[i]]=PGS.TRI(pgs_offspring = PRS_fam[,1], pgs_mother = PRS_fam[,2], pgs_father=PRS_fam[,3], side=2,GxE_int = FALSE)$res_beta
  fit <- glm(D ~ G_logistic, family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_fam_size1[[i]]=summary(fit)$coef
  
  pc_sim = pc_child[c(id1,id2),1:10]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc10 <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_pc10_fam_size1[[i]]=summary(fit_pc10)$coef
  
  
  dat_logit = data.frame(cbind(G_logistic,pc_sim,pheno_sim))
  fit_pc10_geo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_pc10_geo_fam_size1[[i]]=summary(fit_pc10_geo)$coef
  
  
  pc_sim = pheno_mother[c(id1,id2),11:20]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_mo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_mother[[i]]=summary(fit_pc_mo)$coef
  
  pc_sim = pheno_father[c(id1,id2),11:20]
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_fa <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_father[[i]]=summary(fit_pc_fa)$coef
  
  pc_sim = cbind(pheno_father[c(id1,id2),11:20],pheno_mother[c(id1,id2),11:20])
  dat_logit = data.frame(cbind(G_logistic,pc_sim))
  fit_pc_famo <- glm(D ~ ., data = dat_logit,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_pc_famo[[i]]=summary(fit_pc_famo)$coef
  

  
  res_ptdt_size1[[i]]=ptdt(PRS_fam[,1],PRS_fam[,2], PRS_fam[,3],side0=2)$res_beta
  
  
}

save(res_fam_size1,res_glm_fam_size1,res_ptdt_size1,res_pc_father,res_pc_mother,
     res_glm_pc10_fam_size1,res_glm_pc10_geo_fam_size1,res_pc_famo,file=paste0(path,"res_beta04_size3k",".rda"))


