#simulation work for case-parent trio
#Prospectively simulate the data
#GxE investigation
#May 25, 2024
#---------------------------------------------

#Scenario 1: cor_strat=cor_strat_e=0, rho2 (pop strat for GxE term in intercept) =0, no correlations between PRS and E
#beta GE = 0
path<-"/users/zwang4/family/simulation_gxe/05232024/power/e2/"
source("/users/zwang4/family/R/simulation.R")
source("/users/zwang4/family/R/PGS-TRI.R")
# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
task_id <- as.numeric(slurm_arrayid)

print(task_id)
set.seed(task_id)



#Start simulation

#Scenario: type 1 error, beta_prs=0.4, cor_strat=0.25
beta_prs=0.4
betaE_bin_tmp=0
betaE_norm_tmp=0.4
betaGE_normPRS_bin_tmp=0
betaGE_normPRS_norm_tmp=0
fam_size1=200
fam_size2=500
fam_size3=1000
fam_size4=2000
cor_strat1=0
cor_strat2=0

#This part does not change
sim_N=50

res_fam_size1=list()
res_glm_fam_size1=list()
res_caseonly_size1=list()
res_fam_size2=list()
res_glm_fam_size2=list()
res_caseonly_size2=list()
res_fam_size3=list()
res_glm_fam_size3=list()
res_caseonly_size3=list()
res_fam_size4=list()
res_glm_fam_size4=list()
res_caseonly_size4=list()

for(i in 1:sim_N){
  dat=sim_prospective_population_gxe(n_fam=800000,cor_e_prs=F,cor_strat=cor_strat1,cor_strat_e = 0,rho2=cor_strat2,alpha_fam=-6,betaG_normPRS=beta_prs,betaE_bin=betaE_bin_tmp, betaE_norm=betaE_norm_tmp,betaGE_normPRS_bin=betaGE_normPRS_bin_tmp, betaGE_normPRS_norm=betaGE_normPRS_norm_tmp,envir=T)
  
  #First calculate family size 1 = 200
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size1,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size1,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(dat$pgs_fam[id1,1])),as.numeric(as.vector(dat$pgs_fam[id2,1])))
  PRS_fam=dat$pgs_fam[id2,]
  E=rbind(dat$E_sim[id1,],dat$E_sim[id2,])
  
  res_fam_size1[[i]]=PGS.TRI(PRS_fam[,1], PRS_fam[,2], PRS_fam[,3],formula= ~ E_sim_norm,E = dat$E_sim[id2,], GxE_int = T)
  fit <- glm(D ~ G_logistic+E[,2]+G_logistic:E[,2], family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_fam_size1[[i]]=summary(fit)$coef
  
  res_caseonly_size1[[i]]=function_caseonly(prs_c = PRS_fam[,1],envir = dat$E_sim[id2,2],mean_prs=mean(dat$pgs_fam[,1]))
  
  #calculate family size 2 = 500
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size2,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size2,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(dat$pgs_fam[id1,1])),as.numeric(as.vector(dat$pgs_fam[id2,1])))
  PRS_fam=dat$pgs_fam[id2,]
  E=rbind(dat$E_sim[id1,],dat$E_sim[id2,])
  
  res_fam_size2[[i]]=PGS.TRI(PRS_fam[,1], PRS_fam[,2], PRS_fam[,3],formula= ~ E_sim_norm,E = dat$E_sim[id2,], GxE_int = T)
  fit <- glm(D ~ G_logistic+E[,2]+G_logistic:E[,2], family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_fam_size2[[i]]=summary(fit)$coef
  res_caseonly_size2[[i]]=function_caseonly(prs_c = PRS_fam[,1],envir = dat$E_sim[id2,2],mean_prs=mean(dat$pgs_fam[,1]))
  
  #calculate family size 3 = 1000
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size3,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size3,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(dat$pgs_fam[id1,1])),as.numeric(as.vector(dat$pgs_fam[id2,1])))
  PRS_fam=dat$pgs_fam[id2,]
  
  E=rbind(dat$E_sim[id1,],dat$E_sim[id2,])
  
  res_fam_size3[[i]]=PGS.TRI(PRS_fam[,1], PRS_fam[,2], PRS_fam[,3],formula= ~ E_sim_norm,E = dat$E_sim[id2,], GxE_int = T)
  fit <- glm(D ~ G_logistic+E[,2]+G_logistic:E[,2], family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_fam_size3[[i]]=summary(fit)$coef
  res_caseonly_size3[[i]]=function_caseonly(prs_c = PRS_fam[,1],envir = dat$E_sim[id2,2],mean_prs=mean(dat$pgs_fam[,1]))
  
  #calculate family size 4 = 2000
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size4,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size4,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(dat$pgs_fam[id1,1])),as.numeric(as.vector(dat$pgs_fam[id2,1])))
  PRS_fam=dat$pgs_fam[id2,]
  E=rbind(dat$E_sim[id1,],dat$E_sim[id2,])
  
  res_fam_size4[[i]]=PGS.TRI(PRS_fam[,1], PRS_fam[,2], PRS_fam[,3],formula= ~ E_sim_norm,E = dat$E_sim[id2,], GxE_int = T)
  fit <- glm(D ~ G_logistic+E[,2]+G_logistic:E[,2], family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_fam_size4[[i]]=summary(fit)$coef
  res_caseonly_size4[[i]]=function_caseonly(prs_c = PRS_fam[,1],envir = dat$E_sim[id2,2],mean_prs=mean(dat$pgs_fam[,1]))
  
}

save(res_fam_size1,res_glm_fam_size1,res_caseonly_size1,res_fam_size2,res_glm_fam_size2,res_caseonly_size2,
     res_fam_size3,res_glm_fam_size3,res_caseonly_size3,res_fam_size4,res_glm_fam_size4,res_caseonly_size4,file=paste0(path,"power_beta0_",task_id,".rda"))