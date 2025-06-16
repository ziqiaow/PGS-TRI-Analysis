#scenario 3
path<-"/dcs04/nilanjan/data/zwang/family/simulation/snipar/02172025/"
source("/dcs04/nilanjan/data/zwang/family/simulation/snipar/sim_disease_risk_snipar.R")
source("/users/zwang4/family/R/PGS-TRI.R")
source("/users/zwang4/family/R/pTDT.R")
load("/dcs04/nilanjan/data/zwang/family/simulation/snipar/simulated_data/PGS_family.RData")

#Start simulation

#Scenario: type 1 error, beta_prs=0, cor_strat=0
beta_prs=0
cor_strat1=0.25
beta_mother = 0.4
beta_father = 0

fam_size1=1000

PRS_fam = cbind(pgs_children_kid1$SCORE1_SUM,pgs_m$SCORE1_SUM,pgs_f$SCORE1_SUM)

#This part does not change
sim_N=1000

res_fam_size1=list()
res_fam_parent_size1=list()
res_glm_fam_size1=list()
res_ptdt_size1=list()



for(i in 1:sim_N){
  set.seed(i)
  dat=sim_prospective_population(pgs_c = pgs_children_kid1$SCORE1_SUM,pgs_m = pgs_m$SCORE1_SUM,pgs_f = pgs_f$SCORE1_SUM,cor_strat=cor_strat1,alpha_fam=-5,betaG_normPRS=beta_prs,beta_M=beta_mother,beta_F=beta_father)
  
  #First calculate family size 1 = 200
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size1,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size1,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(PRS_fam[id1,1])),as.numeric(as.vector(PRS_fam[id2,1])))
  G_m = c(as.numeric(as.vector(PRS_fam[id1,2])),as.numeric(as.vector(PRS_fam[id2,2])))
  #E = c(as.numeric(as.vector(dat$E_sim[id1])),as.numeric(as.vector(dat$E_sim[id2])))
  PRS_fam_select=PRS_fam[id2,]
  
  
  res_fam_size1[[i]]=PGS.TRI(pgs_offspring = PRS_fam_select[,1], pgs_mother = PRS_fam_select[,2], pgs_father=PRS_fam_select[,3], side=2,GxE_int = FALSE)
  res_fam_parent_size1[[i]]=PGS.TRI(pgs_offspring = PRS_fam_select[,1], pgs_mother = PRS_fam_select[,2], pgs_father=PRS_fam_select[,3], side=2,GxE_int = FALSE, parental_indirect = TRUE)
  
  fit <- glm(D ~ G_logistic+G_m , family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_fam_size1[[i]]=summary(fit)$coef
  res_ptdt_size1[[i]]=ptdt(PRS_fam_select[,1],PRS_fam_select[,2], PRS_fam_select[,3],side0=2)$res_beta
  
}

save(res_fam_size1,res_fam_parent_size1,res_glm_fam_size1,res_ptdt_size1,
     file=paste0(path,"s3_direct0_indirect04",".rda"))




