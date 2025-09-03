#summary of nurturing effect
#Aug 19, 2025
#----------------------------------------------------------

setwd("~/work/family/simulation/snipar/08192025/")
path<-"~/work/family/simulation/snipar/08192025/"


files <- list.files(pattern = "^s.*\\.rda$")
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
  
  res_estimate_fam_size1=array(0,c(sim_N,n_envir))
  res_sd_fam_size1=array(0,c(sim_N,n_envir))
  res_p_fam_size1=array(0,c(sim_N,n_envir))

  res_estimate_glm_fam_size1=array(0,c(sim_N,n_envir))
  res_sd_glm_fam_size1=array(0,c(sim_N,n_envir))
  res_p_glm_fam_size1=array(0,c(sim_N,n_envir))
  

  
  for(j in 1:sim_N){
    res_estimate_fam_size1[j,]=res_fam_parent_size1[[j]]$Coefficients_indirect[,1]
    res_sd_fam_size1[j,]=res_fam_parent_size1[[j]]$Coefficients_indirect[,2]
    res_p_fam_size1[j,]=res_fam_parent_size1[[j]]$Coefficients_indirect[,4]
    

    res_estimate_glm_fam_size1[j,]=res_glm_fam_parent_diff[[j]][,1]
    res_sd_glm_fam_size1[j,]=res_glm_fam_parent_diff[[j]][,2]
    res_p_glm_fam_size1[j,]=res_glm_fam_parent_diff[[j]][,4]
    
    
  }
  tmp=array(0,c(2,n_envir))
  tmp[1,]=apply(res_estimate_fam_size1,2,mean)-true_delta[h]
  tmp[2,]=apply(res_estimate_glm_fam_size1,2,mean)-true_delta[h]
 
  tmp=data.frame(tmp)
  tmp$method=c("fam","logistic")
  colnames(tmp)[1]=c("bias_delta_prs")
  
  
  
  summary_bias[[h]]=tmp
  
  tmp=array(0,c(2,n_envir))
  tmp[1,]=apply(res_sd_fam_size1,2,mean)
  tmp[2,]=apply(res_sd_glm_fam_size1,2,mean)
 
  
  tmp=data.frame(tmp)
  tmp$method=c("fam","logistic")
  colnames(tmp)[1]=c("sd_delta_prs")
  
  
  summary_sd[[h]]=tmp
  
  tmp=array(0,c(2,n_envir))
  tmp[1,]=apply(res_p_fam_size1,2,function(x) mean(x<0.05))
  tmp[2,]=apply(res_p_glm_fam_size1,2,function(x) mean(x<0.05))

  tmp=data.frame(tmp)
  tmp$method=c("fam","logistic")
  colnames(tmp)[1]=c("type1error_delta_prs")
  
  summary_type1error_power[[h]]=tmp
  
  tmp=array(0,c(2,n_envir))
  tmp[1,]=apply(res_estimate_fam_size1,2,sd)
  tmp[2,]=apply(res_estimate_glm_fam_size1,2,sd)
 
  
  tmp=data.frame(tmp)
  tmp$method=c("fam","logistic")
  colnames(tmp)[1]=c("empirical_sd_delta_prs")
  
  summary_empirical_sd[[h]]=tmp
  
  
  
  tmp=array(0,c(2,n_envir))
  tmp[1,]=mean((res_estimate_fam_size1-true_delta[h])^2)
  tmp[2,]=mean((res_estimate_glm_fam_size1-true_delta[h])^2)
  
  
  tmp=data.frame(tmp)
  tmp$method=c("fam","logistic")
  colnames(tmp)[1]=c("mse_delta_prs")
  
  summary_mse[[h]]=tmp
  
  
  #95% coverage probability
  
  tmp=array(0,c(2,n_envir))
  
  
  upper=res_estimate_fam_size1+qnorm(0.975)*res_sd_fam_size1
  lower=res_estimate_fam_size1-qnorm(0.975)*res_sd_fam_size1
  tmp[1,1]=mean(ifelse(true_delta[h]<=upper[,1] & true_delta[h]>=lower[,1],1,0))
  
  
  upper=res_estimate_glm_fam_size1+qnorm(0.975)*res_sd_glm_fam_size1
  lower=res_estimate_glm_fam_size1-qnorm(0.975)*res_sd_glm_fam_size1
  tmp[2,1]=mean(ifelse(true_delta[h]<=upper[,1] & true_delta[h]>=lower[,1],1,0))
  

  tmp=data.frame(tmp)
  tmp$method=c("fam","logistic")
  colnames(tmp)[1]=c("cp_delta_prs")
  
  summary_ci[[h]]=tmp
  
  
}
summary_delta=list()
tmp=data.frame(array(0,c(2,7)))
for(i in 1:4){
  tmp[,1]=summary_bias[[i]][,1]
  tmp[,2]=summary_sd[[i]][,1]
  tmp[,3]=summary_empirical_sd[[i]][,1]
  tmp[,4]=summary_type1error_power[[i]][,1]
  tmp[,5]=summary_ci[[i]][,1]
  tmp[,6]=summary_mse[[i]][,1]
  tmp[,7]=summary_bias[[i]][,2]

  colnames(tmp)=c("bias","sd","empirical_sd","type1error/power","coverage_prob","mse","method")
  summary_delta[[i]]=data.frame(tmp)
  
}

summary_delta[[1]]$true_beta=0
summary_delta[[2]]$true_beta=0
summary_delta[[3]]$true_beta=0.4
summary_delta[[4]]$true_beta=0.4


summary_delta[[1]]$true_delta=0
summary_delta[[2]]$true_delta=0
summary_delta[[3]]$true_delta=0
summary_delta[[4]]$true_delta=0

summary_delta[[1]]$triosize=1000
summary_delta[[2]]$triosize=1000
summary_delta[[3]]$triosize=1000
summary_delta[[4]]$triosize=1000

summary_delta[[1]]$cor=0.4
summary_delta[[2]]$cor=0
summary_delta[[3]]$cor=0.4
summary_delta[[4]]$cor=0


save(summary_delta,file=paste0(path,"results_summary_delta_indirect.rda"))



