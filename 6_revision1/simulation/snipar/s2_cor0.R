#scenario 1
path<-"/dcs04/nilanjan/data/zwang/family/simulation/snipar/08192025/"
source("/dcs04/nilanjan/data/zwang/family/simulation/snipar/sim_disease_risk_snipar.R")
source("/users/zwang4/family/R/PGS-TRI.R")
source("/users/zwang4/family/R/pTDT.R")
load("/dcs04/nilanjan/data/zwang/family/simulation/snipar/simulated_data/PGS_family.RData")
load("/dcs04/nilanjan/data/zwang/family/simulation/snipar/simulated_data/PGS_parents_residuals.RData")
res_mf = 0.5*(res1+res2)

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

#Scenario: beta_prs=0.4, cor_strat=0
beta_prs=0.4
cor_strat1=0
beta_mother = 0
beta_father = 0

fam_size1=1000

PRS_fam = cbind(pgs_children_kid1$SCORE1_SUM,pgs_m$SCORE1_SUM,pgs_f$SCORE1_SUM)

#This part does not change
sim_N=1000

res_fam_size1=res_glm_fam_parent=res_glm_fam_parent_diff=list()
res_fam_parent_size1=list()
res_glm_fam_size1=list()
res_ptdt_size1=list()



for(i in 1:sim_N){
  set.seed(i)
  dat=sim_prospective_population_pheno(pgs_c = pgs_children_kid1$SCORE1_SUM,pgs_m = pgs_m$SCORE1_SUM,pgs_f = pgs_f$SCORE1_SUM,cor_strat=cor_strat1,alpha_fam=-4.9,pheno=res_mf,betaG_normPRS=beta_prs,beta_M=beta_mother,beta_F=beta_father)
  
  
  #First calculate family size 1 = 200
  # randomly select the ncontrol controls and ncase cases
  id1=sample(which(dat$D_sim==0),size=fam_size1,replace=F)
  id2=sample(which(dat$D_sim==1),size=fam_size1,replace=F)
  D = c(dat$D_sim[id1],dat$D_sim[id2])  
  G_logistic = c(as.numeric(as.vector(PRS_fam[id1,1])),as.numeric(as.vector(PRS_fam[id2,1])))
  G_mother = c(as.numeric(as.vector(PRS_fam[id1,2])),as.numeric(as.vector(PRS_fam[id2,2])))
  G_father = c(as.numeric(as.vector(PRS_fam[id1,3])),as.numeric(as.vector(PRS_fam[id2,3])))
  
  #E = c(as.numeric(as.vector(dat$E_sim[id1])),as.numeric(as.vector(dat$E_sim[id2])))
  PRS_fam_select=PRS_fam[id2,]
  
  
  res_fam_size1[[i]]=PGS.TRI(pgs_offspring = PRS_fam_select[,1], pgs_mother = PRS_fam_select[,2], pgs_father=PRS_fam_select[,3], side=2,GxE_int = FALSE)
  res_fam_parent_size1[[i]]=PGS.TRI(pgs_offspring = PRS_fam_select[,1], pgs_mother = PRS_fam_select[,2], pgs_father=PRS_fam_select[,3], side=2,GxE_int = FALSE, parental_indirect = TRUE)
  
  fit <- glm(D ~ G_logistic , family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_fam_size1[[i]]=summary(fit)$coef
  
  fit2 <- glm(D ~ G_logistic+G_mother+G_father, family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm_fam_parent[[i]] = summary(fit2)$coef
  res_glm_fam_parent_diff[[i]] = difftest_glm("G_mother","G_father",fit2)
  
  res_ptdt_size1[[i]]=ptdt(PRS_fam_select[,1],PRS_fam_select[,2], PRS_fam_select[,3],side0=2)$statistic
  
  
}

save(res_fam_size1,res_fam_parent_size1,res_glm_fam_size1,res_ptdt_size1,res_glm_fam_parent,res_glm_fam_parent_diff,
     file=paste0(path,"s2_direct04_indirect0_cor0",".rda"))

