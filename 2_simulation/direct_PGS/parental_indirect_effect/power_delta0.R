#Example code for power analyses in the manuscript
#When delta-IDE = 0

path<-"/users/zwang4/family/simulation/nurture/redo/05212024/power/"
source("/users/zwang4/family/R/simulation.R")
source("/users/zwang4/family/R/PGS-TRI.R")

# grab the array id value from the environment variable passed from sbatch
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')

# coerce the value to an integer
task_id <- as.numeric(slurm_arrayid)

print(task_id)
set.seed(task_id)

#write my own code
difftest_glm <- function(x1, x2, model){
  diffest <- summary(model)$coef[x1,"Estimate"]-summary(model)$coef[x2,"Estimate"]
  vardiff <- (summary(model)$coef[x1,"Std. Error"]^2 + 
                summary(model)$coef[x2,"Std. Error"]^2) - (2*(vcov(model)[x1, x2])) 
  # variance of x1 + variance of x2 - 2*covariance of x1 and x2
  diffse <- sqrt(vardiff)
  zdiff <- (diffest)/(diffse)
  ptdiff <- 2* pnorm(abs(zdiff), lower.tail = FALSE)
  upr <- diffest + qnorm(0.975)*diffse # will usually be very close to 1.96
  lwr <- diffest - qnorm(0.975)*diffse
  
  res_delta=array(0,c(1,6))
  res_delta[1,]=c(diffest,diffse,zdiff,ptdiff,lwr,upr)
  colnames(res_delta)=c("Estimate","Std.Error","Z.value","Pvalue","lwr","upr")
  rownames(res_delta)="Nurture_Diff_MF"
  return(res_delta)
}

#Start simulation

#Scenario: type 1 error, beta_prs=0, cor_strat=0
beta_prs=0
cor_strat1=0
beta_M_tmp=0
beta_F_tmp=0

fam_size1=200
fam_size2=500
fam_size3=1000
fam_size4=2000


#This part does not change
sim_N=50

res_fam_size1=list()
res_fam_size2=list()
res_fam_size3=list()
res_fam_size4=list()


res_glm_fam_size1=list()
res_glm_fam_size2=list()
res_glm_fam_size3=list()
res_glm_fam_size4=list()


for(i in 1:sim_N){
  dat=sim_prospective_population(n_fam=500000,cor_strat=cor_strat1,alpha_fam=-5.5,betaG_normPRS=beta_prs,envir=FALSE,nurture=TRUE,beta_M=beta_M_tmp,beta_F=beta_F_tmp)
  
  #First calculate family size 1 = 200
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size1,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size1,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(dat$pgs_fam[id1,1])),as.numeric(as.vector(dat$pgs_fam[id2,1])))
  G_mother = c(as.numeric(as.vector(dat$pgs_fam[id1,2])),as.numeric(as.vector(dat$pgs_fam[id2,2])))
  G_father = c(as.numeric(as.vector(dat$pgs_fam[id1,3])),as.numeric(as.vector(dat$pgs_fam[id2,3])))
  
  PRS_fam=dat$pgs_fam[id2,]
  
  res_fam_size1[[i]]=PGS.TRI(pgs_offspring = PRS_fam[,1], pgs_mother = PRS_fam[,2], pgs_father=PRS_fam[,3], side=2,GxE_int = FALSE,parental_indirect =T)
  fit=  glm(D ~ G_logistic+G_mother+G_father, family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  
  glm_tmp=list()
  glm_tmp$coef=summary(fit)$coef
  glm_tmp$res_delta=difftest_glm("G_mother","G_father",fit)
  res_glm_fam_size1[[i]]=glm_tmp
  
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size2,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size2,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(dat$pgs_fam[id1,1])),as.numeric(as.vector(dat$pgs_fam[id2,1])))
  PRS_fam=dat$pgs_fam[id2,]
  G_mother = c(as.numeric(as.vector(dat$pgs_fam[id1,2])),as.numeric(as.vector(dat$pgs_fam[id2,2])))
  G_father = c(as.numeric(as.vector(dat$pgs_fam[id1,3])),as.numeric(as.vector(dat$pgs_fam[id2,3])))
  
  res_fam_size2[[i]]=PGS.TRI(pgs_offspring = PRS_fam[,1], pgs_mother = PRS_fam[,2], pgs_father=PRS_fam[,3], side=2,GxE_int = FALSE,parental_indirect =T)
  fit=  glm(D ~ G_logistic+G_mother+G_father, family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  
  glm_tmp=list()
  glm_tmp$coef=summary(fit)$coef
  glm_tmp$res_delta=difftest_glm("G_mother","G_father",fit)
  res_glm_fam_size2[[i]]=glm_tmp
  
  #calculate family size 3 = 1000
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size3,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size3,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(dat$pgs_fam[id1,1])),as.numeric(as.vector(dat$pgs_fam[id2,1])))
  PRS_fam=dat$pgs_fam[id2,]
  G_mother = c(as.numeric(as.vector(dat$pgs_fam[id1,2])),as.numeric(as.vector(dat$pgs_fam[id2,2])))
  G_father = c(as.numeric(as.vector(dat$pgs_fam[id1,3])),as.numeric(as.vector(dat$pgs_fam[id2,3])))
  res_fam_size3[[i]]=PGS.TRI(pgs_offspring = PRS_fam[,1], pgs_mother = PRS_fam[,2], pgs_father=PRS_fam[,3], side=2,GxE_int = FALSE,parental_indirect =T)
  fit=  glm(D ~ G_logistic+G_mother+G_father, family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  
  glm_tmp=list()
  glm_tmp$coef=summary(fit)$coef
  glm_tmp$res_delta=difftest_glm("G_mother","G_father",fit)
  res_glm_fam_size3[[i]]=glm_tmp
  
  
  #calculate family size 4 = 2000
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size4,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size4,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(dat$pgs_fam[id1,1])),as.numeric(as.vector(dat$pgs_fam[id2,1])))
  PRS_fam=dat$pgs_fam[id2,]
  G_mother = c(as.numeric(as.vector(dat$pgs_fam[id1,2])),as.numeric(as.vector(dat$pgs_fam[id2,2])))
  G_father = c(as.numeric(as.vector(dat$pgs_fam[id1,3])),as.numeric(as.vector(dat$pgs_fam[id2,3])))
  res_fam_size4[[i]]=PGS.TRI(pgs_offspring = PRS_fam[,1], pgs_mother = PRS_fam[,2], pgs_father=PRS_fam[,3], side=2,GxE_int = FALSE,parental_indirect =T)
  fit=  glm(D ~ G_logistic+G_mother+G_father, family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  
  glm_tmp=list()
  glm_tmp$coef=summary(fit)$coef
  glm_tmp$res_delta=difftest_glm("G_mother","G_father",fit)
  res_glm_fam_size4[[i]]=glm_tmp
  
  
  
}

save(res_fam_size1,res_glm_fam_size1,
     res_fam_size2,res_glm_fam_size2,
     res_fam_size3,res_glm_fam_size3,
     res_fam_size4,res_glm_fam_size4,file=paste0(path,"power_delta0_",task_id,".rda"))



