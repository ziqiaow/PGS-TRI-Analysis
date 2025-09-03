#Next simulate disease status using population disease risk model
#Data based on 20 generations of assortative mating of 100,000 families
#Program: https://snipar.readthedocs.io/en/latest/simulation.html
#Reference: https://www.nature.com/articles/s41588-022-01085-0
#Family-based PGS-TRI
#Likelihood and parameter estimation
#Feb 17, 2025
#--------------------------------------------------

sim_prospective_population_pheno=function(pgs_c,pgs_m,pgs_f,cor_strat=0.5,alpha_fam=-5,betaG_normPRS=0.4,beta_M=0,beta_F=0,pheno){
  

  alpha_sim = rnorm(length(pgs_c),mean = alpha_fam+cor_strat*pheno, sd = 1)
  
  xb <- alpha_sim + betaG_normPRS*pgs_c + beta_M*pgs_m + beta_F*pgs_f
  
  prob <- 1/(1 + exp(-xb)) #Even though our assumption is that the disease model follows log-linear, we simulate data using logistic regression
  D_sim <- rbinom(n = length(pgs_c), size = 1, prob = prob)  # 1=case, 0=control
  cat("Disease prevalence: ",mean(D_sim))
  
  #### Bind D, G_logistic,PRS_fam, and E into a list
  dat = list(D_sim=D_sim)#, pgs_fam=pgs_fam)
  return(dat)
  
  
}
