#Simulation of selection bias
#--------------------------------------------------
#selection bias in direct effect
sim_prospective_population_select=function(n_fam=100000,mu_select = -1, beta_select_E1=0.4,beta_select_E2=0, cor_e_prs=F,cor_strat=0.25,alpha_fam=-5.5,betaG_normPRS=0.4,betaE1=0.4, betaE2=0.4, PGScauseE = FALSE){
  library(MASS)
  # random generation for sigma_i^2 of the PGS values for the family, assume a mixture of gamma distributions
  rmixgamma <- function(n, pi, alpha, beta) {
    k <- sample.int(length(pi), n, replace = TRUE, prob = pi)
    rgamma(n, alpha[k], beta[k])
  }
  
  pi <- c(0.6,0.3,0.1)
  alpha <- c(6, 15,60)
  beta <- c(15, 30,150)
  
  var_sim=rmixgamma(n_fam, pi, alpha, beta)
  if(PGScauseE == TRUE) {
    
    mu_sim=c(0,alpha_fam) #this is the mean vector for family PRS mu_i and alpha_i in the disease model
    sigma_sim=matrix(c(1,cor_strat,cor_strat,1),2,2) #assume there is population stratification
    mu_alpha=MASS::mvrnorm(n=n_fam,mu=mu_sim,Sigma=sigma_sim)
    mu_var=cbind(mu_alpha[,1],var_sim)
    myfunctionpgs=function(x){
      sigma_tmp=matrix(c(1,0.5,0.5,0.5,1,0,0.5,0,1)*x[2],3,3)
      MASS::mvrnorm(n=1,mu=rep(x[1],3),Sigma=sigma_tmp)
      
    }
    pgs_fam=t(apply(mu_var,1,myfunctionpgs))
    
    colnames(pgs_fam)=c("pgs_c","pgs_m","pgs_f")
    prs <- as.numeric(as.vector(pgs_fam[,1]))
    
    E1 = 0.2*prs + rnorm(n_fam)
    
    xb <- mu_alpha[,2] + betaG_normPRS*prs + betaE1*E1
    
    prob <- 1/(1 + exp(-xb)) #Even though our assumption is that the disease model follows log-linear, we simulate data using logistic regression
    D_sim <- rbinom(n = n_fam, size = 1, prob = prob)  # 1=case, 0=control
    cat("Disease prevalence: ",mean(D_sim))
    
    prob_select <- 1/(1+exp(-beta_select_E1 * E1 - mu_select)) #if mu_select =0, beta_select=0, prob_select = 0.5 (randomly selected)
    select_sim = rbinom(n = n_fam, size =1, prob = prob_select)
    
    
    id = which(select_sim == 1)
    D_sim_select = D_sim[id]
    E_sim_norm_select = E1[id]
    pgs_fam_select = pgs_fam[id,]
    
    cat("Select",length(id),"families")
    #### Bind D, G_logistic,PRS_fam, and E into a list
    dat = list(D_sim=D_sim_select, E_sim=E_sim_norm_select, pgs_fam=pgs_fam_select)
    
    
    
  } else {
    if(cor_e_prs == TRUE){
      
      rho1=0.2
      rho2=0
      mu1=rnorm(n_fam)
      e1=rnorm(n_fam)
      e2=rnorm(n_fam)
      mu_alpha=array(0,c(n_fam,3))
      for(i in 1:n_fam){
        sigma_sim=matrix(c(1,rho1,rho2,rho1,1,0,rho2,0,1),3,3)
        sigma_sim=chol(sigma_sim)
        mu_alpha[i,]=c(mu1[i],e1[i],e2[i]) %*% t(sigma_sim) 
      }
      
      
      mu_var=cbind(mu_alpha[,1],var_sim)
      myfunctionpgs=function(x){
        sigma_tmp=matrix(c(1,0.5,0.5,0.5,1,0,0.5,0,1)*x[2],3,3)
        MASS::mvrnorm(n=1,mu=rep(x[1],3),Sigma=sigma_tmp)
        
      }
      pgs_fam=t(apply(mu_var,1,myfunctionpgs))
      
      colnames(pgs_fam)=c("pgs_c","pgs_m","pgs_f")
      prs <- as.numeric(as.vector(pgs_fam[,1]))
      
      E1 <- mu_alpha[,2]
      E2 <- mu_alpha[,3]
      
      
      strat=rnorm(n=n_fam,mean=alpha_fam+cor_strat*mu_alpha[,1], sd=1)
      
      xb <- strat + betaG_normPRS*prs + betaE1*E1  + betaE2*E2
     
      prob <- 1/(1 + exp(-xb)) #Even though our assumption is that the disease model follows log-linear, we simulate data using logistic regression
      D_sim <- rbinom(n = n_fam, size = 1, prob = prob)  # 1=case, 0=control
      cat("Disease prevalence: ",mean(D_sim))
      
      prob_select <- 1/(1+exp(-beta_select_E1 * E1 -beta_select_E2 * E2 - mu_select)) #if mu_select =0, beta_select=0, prob_select = 0.5 (randomly selected)
      select_sim = rbinom(n = n_fam, size =1, prob = prob_select)
      
      id = which(select_sim == 1)
      D_sim_select = D_sim[id]
   
      pgs_fam_select = pgs_fam[id,]
      
      if(beta_select_E2 == 0){
        E_sim_norm_select = E1[id]
      } else {
        E_sim_norm_select = cbind(E1[id],E2[id])
      }
     
      cat("Select",length(id),"families")
      #### Bind D, G_logistic,PRS_fam, and E into a list
      dat = list(D_sim=D_sim_select, E_sim=E_sim_norm_select, pgs_fam=pgs_fam_select)
      
    }else{
      
    E1=rnorm(n_fam)
    mu_sim=c(0,alpha_fam) #this is the mean vector for family PRS mu_i and alpha_i in the disease model
    sigma_sim=matrix(c(1,cor_strat,cor_strat,1),2,2) #assume there is population stratification
    mu_alpha=MASS::mvrnorm(n=n_fam,mu=mu_sim,Sigma=sigma_sim)
    mu_var=cbind(mu_alpha[,1],var_sim)
    myfunctionpgs=function(x){
      sigma_tmp=matrix(c(1,0.5,0.5,0.5,1,0,0.5,0,1)*x[2],3,3)
      MASS::mvrnorm(n=1,mu=rep(x[1],3),Sigma=sigma_tmp)
      
    }
    pgs_fam=t(apply(mu_var,1,myfunctionpgs))
    
    colnames(pgs_fam)=c("pgs_c","pgs_m","pgs_f")
    prs <- as.numeric(as.vector(pgs_fam[,1]))
    
    xb <- mu_alpha[,2] + betaG_normPRS*prs + betaE1*E1
    
    prob <- 1/(1 + exp(-xb)) #Even though our assumption is that the disease model follows log-linear, we simulate data using logistic regression
    D_sim <- rbinom(n = n_fam, size = 1, prob = prob)  # 1=case, 0=control
    cat("Disease prevalence: ",mean(D_sim))
    
    prob_select <- 1/(1+exp(-beta_select_E1 * E1 - mu_select)) #if mu_select =0, beta_select=0, prob_select = 0.5 (randomly selected)
    select_sim = rbinom(n = n_fam, size =1, prob = prob_select)
    
    
    id = which(select_sim == 1)
    D_sim_select = D_sim[id]
    E_sim_norm_select = E1[id]
    pgs_fam_select = pgs_fam[id,]
    
    cat("Select",length(id),"families")
    #### Bind D, G_logistic,PRS_fam, and E into a list
    dat = list(D_sim=D_sim_select, E_sim=E_sim_norm_select, pgs_fam=pgs_fam_select)
    
    
    }
    
    
  }
  
  return(dat)
  
  
}






#selection bias in direct effect and parental indirect effect
sim_prospective_population_parents=function(n_fam=100000,mu_select = -1,beta_select_m = 0.2 ,beta_select_f = 0.2 , beta_select_E1=0.4, cor_strat=0.25,alpha_fam=-5.5,betaG_normPRS=0.4,betaE1=0.4 , envir = FALSE){
  library(MASS)
  # random generation for sigma_i^2 of the PGS values for the family, assume a mixture of gamma distributions
  rmixgamma <- function(n, pi, alpha, beta) {
    k <- sample.int(length(pi), n, replace = TRUE, prob = pi)
    rgamma(n, alpha[k], beta[k])
  }
  
  pi <- c(0.6,0.3,0.1)
  alpha <- c(6, 15,60)
  beta <- c(15, 30,150)
  
  var_sim=rmixgamma(n_fam, pi, alpha, beta)
  if(envir == TRUE) {
    
    mu_sim=c(0,alpha_fam) #this is the mean vector for family PRS mu_i and alpha_i in the disease model
    sigma_sim=matrix(c(1,cor_strat,cor_strat,1),2,2) #assume there is population stratification
    mu_alpha=MASS::mvrnorm(n=n_fam,mu=mu_sim,Sigma=sigma_sim)
    mu_var=cbind(mu_alpha[,1],var_sim)
    myfunctionpgs=function(x){
      sigma_tmp=matrix(c(1,0.5,0.5,0.5,1,0,0.5,0,1)*x[2],3,3)
      MASS::mvrnorm(n=1,mu=rep(x[1],3),Sigma=sigma_tmp)
      
    }
    pgs_fam=t(apply(mu_var,1,myfunctionpgs))
    
    colnames(pgs_fam)=c("pgs_c","pgs_m","pgs_f")
    prs <- as.numeric(as.vector(pgs_fam[,1]))
    
    E1 = rnorm(n_fam)
    
    xb <- mu_alpha[,2] + betaG_normPRS*prs + betaE1*E1
    
    prob <- 1/(1 + exp(-xb)) #Even though our assumption is that the disease model follows log-linear, we simulate data using logistic regression
    D_sim <- rbinom(n = n_fam, size = 1, prob = prob)  # 1=case, 0=control
    cat("Disease prevalence: ",mean(D_sim))
    
    prob_select <- 1/(1+exp(-beta_select_E1 * E1 -beta_select_m * pgs_fam[,2] -beta_select_f * pgs_fam[,3] - mu_select)) #if mu_select =0, beta_select=0, prob_select = 0.5 (randomly selected)
    select_sim = rbinom(n = n_fam, size =1, prob = prob_select)
    
    
    id = which(select_sim == 1)
    D_sim_select = D_sim[id]
    E_sim_norm_select = E1[id]
    pgs_fam_select = pgs_fam[id,]
    
    cat("Select",length(id),"families")
    #### Bind D, G_logistic,PRS_fam, and E into a list
    dat = list(D_sim=D_sim_select, E_sim=E_sim_norm_select, pgs_fam=pgs_fam_select)
    
    
    
  } else {
    
    
    mu_sim=c(0,alpha_fam) #this is the mean vector for family PRS mu_i and alpha_i in the disease model
    sigma_sim=matrix(c(1,cor_strat,cor_strat,1),2,2) #assume there is population stratification
    mu_alpha=MASS::mvrnorm(n=n_fam,mu=mu_sim,Sigma=sigma_sim)
    mu_var=cbind(mu_alpha[,1],var_sim)
    myfunctionpgs=function(x){
      sigma_tmp=matrix(c(1,0.5,0.5,0.5,1,0,0.5,0,1)*x[2],3,3)
      MASS::mvrnorm(n=1,mu=rep(x[1],3),Sigma=sigma_tmp)
      
    }
    pgs_fam=t(apply(mu_var,1,myfunctionpgs))
    
    colnames(pgs_fam)=c("pgs_c","pgs_m","pgs_f")
    prs <- as.numeric(as.vector(pgs_fam[,1]))
    
    
    xb <- mu_alpha[,2] + betaG_normPRS*prs 
    
    prob <- 1/(1 + exp(-xb)) #Even though our assumption is that the disease model follows log-linear, we simulate data using logistic regression
    D_sim <- rbinom(n = n_fam, size = 1, prob = prob)  # 1=case, 0=control
    cat("Disease prevalence: ",mean(D_sim))
    
    prob_select <- 1/(1+exp(-beta_select_m * pgs_fam[,2] -beta_select_f * pgs_fam[,3] - mu_select)) #if mu_select =0, beta_select=0, prob_select = 0.5 (randomly selected)
    select_sim = rbinom(n = n_fam, size =1, prob = prob_select)
    
    
    id = which(select_sim == 1)
    D_sim_select = D_sim[id]
    pgs_fam_select = pgs_fam[id,]
    
    cat("Select",length(id),"families")
    #### Bind D, G_logistic,PRS_fam, and E into a list
    dat = list(D_sim=D_sim_select, pgs_fam=pgs_fam_select)
    
    
    
    
  }
  
  return(dat)
  
  
}



