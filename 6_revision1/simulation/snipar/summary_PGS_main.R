
setwd("~/work/family/simulation/snipar/08192025/")
path<-"~/work/family/simulation/snipar/08192025/"


files <- list.files(pattern = "^s.*\\.rda$")
true_beta=c(0,0,0.4,0.4)
true_delta = c(0,0,0,0)


sim_N=1000
n_envir=1
summary_bias=list()
summary_sd=list()
summary_empirical_sd=list()
summary_type1error_power=list()
summary_mse=list()
summary_ci=list()


h=0
for(i in files){
  
  load(i)
  h=h+1 
  
  #summarize
  n_envir = 1
  
  res_estimate_fam_size1=array(0,c(sim_N,n_envir))
  res_sd_fam_size1=array(0,c(sim_N,n_envir))
  res_p_fam_size1=array(0,c(sim_N,n_envir))
  
  res_estimate_fam_parent_size1=array(0,c(sim_N,n_envir))
  res_sd_fam_parent_size1=array(0,c(sim_N,n_envir))
  res_p_fam_parent_size1=array(0,c(sim_N,n_envir))
  
  res_estimate_glm_fam_size1=array(0,c(sim_N,n_envir))
  res_sd_glm_fam_size1=array(0,c(sim_N,n_envir))
  res_p_glm_fam_size1=array(0,c(sim_N,n_envir))
  
  res_estimate_glm_fam_parent_size1=array(0,c(sim_N,n_envir))
  res_sd_glm_fam_parent_size1=array(0,c(sim_N,n_envir))
  res_p_glm_fam_parent_size1=array(0,c(sim_N,n_envir))
  
  res_estimate_ptdt_fam_size1=array(0,c(sim_N,n_envir))
  res_sd_ptdt_fam_size1=array(0,c(sim_N,n_envir))
  res_p_ptdt_fam_size1=array(0,c(sim_N,n_envir))
  
  
  
  for(j in 1:sim_N){
    res_estimate_fam_size1[j,]=res_fam_size1[[j]]$Coefficients_direct[,1]
    res_sd_fam_size1[j,]=res_fam_size1[[j]]$Coefficients_direct[,2]
    res_p_fam_size1[j,]=res_fam_size1[[j]]$Coefficients_direct[,4]
    
    res_estimate_ptdt_fam_size1[j,]=res_ptdt_size1[[j]][,1]
    res_sd_ptdt_fam_size1[j,]=res_ptdt_size1[[j]][,2]
    res_p_ptdt_fam_size1[j,]=res_ptdt_size1[[j]][,4]
    
    res_estimate_glm_fam_size1[j,]=res_glm_fam_size1[[j]]["G_logistic",1]
    res_sd_glm_fam_size1[j,]=res_glm_fam_size1[[j]]["G_logistic",2]
    res_p_glm_fam_size1[j,]=res_glm_fam_size1[[j]]["G_logistic",4]
    
    res_estimate_glm_fam_parent_size1[j,]=res_glm_fam_parent[[j]]["G_logistic",1]
    res_sd_glm_fam_parent_size1[j,]=res_glm_fam_parent[[j]]["G_logistic",2]
    res_p_glm_fam_parent_size1[j,]=res_glm_fam_parent[[j]]["G_logistic",4]
    
    res_estimate_fam_parent_size1[j,]=res_fam_parent_size1[[j]]$Coefficients_direct[,1]
    res_sd_fam_parent_size1[j,]=res_fam_parent_size1[[j]]$Coefficients_direct[,2]
    res_p_fam_parent_size1[j,]=res_fam_parent_size1[[j]]$Coefficients_direct[,4]
    
    
    
  }
  
  
  tmp=array(0,c(5,n_envir))
  
  
  tmp[1,]=apply(res_estimate_fam_size1,2,mean)-true_beta[h]
  tmp[2,]=apply(res_estimate_glm_fam_size1,2,mean)-true_beta[h]
  tmp[3,]=apply(res_estimate_fam_parent_size1,2,mean)-true_beta[h]
  tmp[4,]=apply(res_estimate_glm_fam_parent_size1,2,mean)-true_beta[h]
  tmp[5,]=apply(res_estimate_ptdt_fam_size1,2,mean)-true_beta[h]
  
  
  tmp=data.frame(tmp)
  tmp$method=c("fam","logistic","fam_parent_adj","logistic_parent_adj","ptdt")
  #tmp$samplesize=c(rep(200,3),rep(500,3),rep(1000,3),rep(2500,3))
  colnames(tmp)[1]=c("bias_beta_prs")
  
  summary_bias[[h]]=tmp
  
  tmp=array(0,c(5,n_envir))
  tmp[1,]=apply(res_sd_fam_size1,2,mean)
  tmp[2,]=apply(res_sd_glm_fam_size1,2,mean)
  tmp[3,]=apply(res_sd_fam_parent_size1,2,mean)
  tmp[4,]=apply(res_sd_glm_fam_parent_size1,2,mean)
  tmp[5,]=apply(res_sd_ptdt_fam_size1,2,mean)
  
  
  tmp=data.frame(tmp)
  tmp$method=c("fam","logistic","fam_parent_adj","logistic_parent_adj","ptdt")
  #tmp$samplesize=c(rep(200,3),rep(500,3),rep(1000,3),rep(2500,3))
  colnames(tmp)[1]=c("sd_beta_prs")
  
  summary_sd[[h]]=tmp
  
  tmp=array(0,c(5,n_envir))
  tmp[1,]=apply(res_p_fam_size1,2,function(x) mean(x<0.05))
  tmp[2,]=apply(res_p_glm_fam_size1,2,function(x) mean(x<0.05))
  tmp[3,]=apply(res_p_fam_parent_size1,2,function(x) mean(x<0.05))
  tmp[4,]=apply(res_p_glm_fam_parent_size1,2,function(x) mean(x<0.05))
  tmp[5,]=apply(res_p_ptdt_fam_size1,2,function(x) mean(x<0.05))
  
  
  
  
  tmp=data.frame(tmp)
  tmp$method=c("fam","logistic","fam_parent_adj","logistic_parent_adj","ptdt")
  colnames(tmp)[1]=c("type1error_beta_prs")
  
  summary_type1error_power[[h]]=tmp
  
  tmp=array(0,c(5,n_envir))
  tmp[1,]=apply(res_estimate_fam_size1,2,sd)
  tmp[2,]=apply(res_estimate_glm_fam_size1,2,sd)
  tmp[3,]=apply(res_estimate_fam_parent_size1,2,sd)
  tmp[4,]=apply(res_estimate_glm_fam_parent_size1,2,sd)
  tmp[5,]=apply(res_estimate_ptdt_fam_size1,2,sd)
  
  
  tmp=data.frame(tmp)
  tmp$method=c("fam","logistic","fam_parent_adj","logistic_parent_adj","ptdt")
  colnames(tmp)[1]=c("empirical_sd_beta_prs")
  
  summary_empirical_sd[[h]]=tmp
  
  tmp=array(0,c(5,n_envir))
  tmp[1,]=mean((res_estimate_fam_size1-true_beta[h])^2)
  tmp[2,]=mean((res_estimate_glm_fam_size1-true_beta[h])^2)
  tmp[3,]=mean((res_estimate_fam_parent_size1-true_beta[h])^2)
  tmp[4,]=mean((res_estimate_glm_fam_parent_size1-true_beta[h])^2)
  tmp[5,]=mean((res_estimate_ptdt_fam_size1-true_beta[h])^2)
  
  
  
  tmp=data.frame(tmp)
  tmp$method=c("fam","logistic","fam_parent_adj","logistic_parent_adj","ptdt")
  colnames(tmp)[1]=c("mse_beta_prs")
  
  summary_mse[[h]]=tmp
  
  
  #95% coverage probability
  
  tmp=array(0,c(5,n_envir))
  
  upper=res_estimate_fam_size1+qnorm(0.975)*res_sd_fam_size1
  lower=res_estimate_fam_size1-qnorm(0.975)*res_sd_fam_size1
  tmp[1,1]=mean(ifelse(true_beta[h]<=upper[,1] & true_beta[h]>=lower[,1],1,0))
  
  
  upper=res_estimate_glm_fam_size1+qnorm(0.975)*res_sd_glm_fam_size1
  lower=res_estimate_glm_fam_size1-qnorm(0.975)*res_sd_glm_fam_size1
  tmp[2,1]=mean(ifelse(true_beta[h]<=upper[,1] & true_beta[h]>=lower[,1],1,0))
  
  
  upper=res_estimate_fam_parent_size1+qnorm(0.975)*res_sd_fam_parent_size1
  lower=res_estimate_fam_parent_size1-qnorm(0.975)*res_sd_fam_parent_size1
  tmp[3,1]=mean(ifelse(true_beta[h]<=upper[,1] & true_beta[h]>=lower[,1],1,0))
  
  
  upper=res_estimate_glm_fam_parent_size1+qnorm(0.975)*res_sd_glm_fam_parent_size1
  lower=res_estimate_glm_fam_parent_size1-qnorm(0.975)*res_sd_glm_fam_parent_size1
  tmp[4,1]=mean(ifelse(true_beta[h]<=upper[,1] & true_beta[h]>=lower[,1],1,0))
  
  upper=res_estimate_ptdt_fam_size1+qnorm(0.975)*res_sd_ptdt_fam_size1
  lower=res_estimate_ptdt_fam_size1-qnorm(0.975)*res_sd_ptdt_fam_size1
  tmp[5,1]=mean(ifelse(true_beta[h]<=upper[,1] & true_beta[h]>=lower[,1],1,0))
  
  
  
  tmp=data.frame(tmp)
  tmp$method=c("fam","logistic","fam_parent_adj","logistic_parent_adj","ptdt")
  colnames(tmp)[1]=c("cp_beta_prs")
  
  summary_ci[[h]]=tmp
}


summary_beta=list()
tmp=data.frame(array(0,c(5,7)))
for(i in 1:4){
  tmp[,1]=summary_bias[[i]][,1]
  tmp[,2]=summary_sd[[i]][,1]
  tmp[,3]=summary_empirical_sd[[i]][,1]
  tmp[,4]=summary_type1error_power[[i]][,1]
  tmp[,5]=summary_ci[[i]][,1]
  tmp[,6]=summary_mse[[i]][,1]
  tmp[,7]=summary_bias[[i]][,2]
  
  colnames(tmp)=c("bias","sd","empirical_sd","type1error/power","coverage_prob","mse","method")
  summary_beta[[i]]=data.frame(tmp)
  
}
summary_beta[[1]]$true_beta=0
summary_beta[[2]]$true_beta=0
summary_beta[[3]]$true_beta=0.4
summary_beta[[4]]$true_beta=0.4


summary_beta[[1]]$triosize=1000
summary_beta[[2]]$triosize=1000
summary_beta[[3]]$triosize=1000
summary_beta[[4]]$triosize=1000

summary_beta[[1]]$cor=0.4
summary_beta[[2]]$cor=0
summary_beta[[3]]$cor=0.4
summary_beta[[4]]$cor=0


save(summary_beta,file=paste0(path,"results_summary_beta_main.rda"))
