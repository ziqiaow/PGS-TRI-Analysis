#Next simulate disease status using population disease risk model
#Family-based PGS-TRI
#Likelihood and parameter estimation
#May 5, 2024
#--------------------------------------------------
sim_prospective_population=function(pgs_c,pgs_m,pgs_f,cor_strat=0.5,alpha_fam=-5,betaG_normPRS=0.4,e_strat,pc=F,mu_pc,envir=F){
    library(MASS)
  
    pgs_fam = cbind(pgs_c,pgs_m,pgs_f)
    colnames(pgs_fam)=c("pgs_c","pgs_m","pgs_f")
    prs <- as.numeric(as.vector(pgs_fam[,1]))
    
    
    if(pc == T){
      alpha_sim = rnorm(length(prs),mean = alpha_fam + cor_strat*e_strat + 0.5*mu_pc, sd=1)
    }else{
    alpha_sim = rnorm(length(prs),mean = alpha_fam + cor_strat*e_strat, sd=1)
    }
    
    xb <- alpha_sim + betaG_normPRS*prs
  
  
  prob <- 1/(1 + exp(-xb)) #Even though our assumption is that the disease model follows log-linear, we simulate data using logistic regression
  D_sim <- rbinom(n = length(prs), size = 1, prob = prob)  # 1=case, 0=control
  cat("Disease prevalence: ",mean(D_sim))
  
  #### Bind D, G_logistic,PRS_fam, and E into a list
  dat = list(D_sim=D_sim, pgs_fam=pgs_fam)
  return(dat)
  
  
}
