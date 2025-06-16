#scenario 8: no population stratification bias
#PGS,E -> D, PGS_m,PGS_f,E-> Z (selection), mother father different selection magnitude
path<-"/dcs04/nilanjan/data/zwang/family/simulation/collider_bias/"
source("/dcs04/nilanjan/data/zwang/family/simulation/collider_bias/simulation_collider.R")
source("/users/zwang4/family/R/PGS-TRI.R")
source("/users/zwang4/family/R/pTDT.R")


#Start simulation

#Scenario: type 1 error, beta_prs=0, cor_strat=0
beta_prs=0.4
cor_strat1=0

fam_size1=1000


#This part does not change
sim_N=1000

res_fam_size1=res_fam_parent_size1=list()
res_glm_fam_size1=list()
res_ptdt_size1=list()




for(i in 1:sim_N){
  set.seed(i)
  dat=sim_prospective_population_parents(n_fam=600000,beta_select_m = 0.2 ,beta_select_f = 0, cor_strat=cor_strat1,alpha_fam=-5.2,betaG_normPRS=beta_prs, beta_select_E1=0.4,envir = T)
  
  #First calculate family size 1 = 200
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size1,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size1,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(dat$pgs_fam[id1,1])),as.numeric(as.vector(dat$pgs_fam[id2,1])))
  #E = c(as.numeric(as.vector(dat$E_sim[id1])),as.numeric(as.vector(dat$E_sim[id2])))
  PRS_fam=dat$pgs_fam[id2,]
  
  
  res_fam_size1[[i]]=PGS.TRI(pgs_offspring = PRS_fam[,1], pgs_mother = PRS_fam[,2], pgs_father=PRS_fam[,3], side=2,GxE_int = FALSE)
  res_fam_parent_size1[[i]]=PGS.TRI(pgs_offspring = PRS_fam[,1], pgs_mother = PRS_fam[,2], pgs_father=PRS_fam[,3], side=2,GxE_int = FALSE, parental_indirect = TRUE)
  
  fit <- glm(D ~ G_logistic , family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_fam_size1[[i]]=summary(fit)$coef
  res_ptdt_size1[[i]]=ptdt(PRS_fam[,1],PRS_fam[,2], PRS_fam[,3],side0=2)$res_beta
  
  
}

save(res_fam_size1,res_glm_fam_size1,res_ptdt_size1,res_fam_parent_size1,
     file=paste0(path,"s8_pgs04_cor0",".rda"))




