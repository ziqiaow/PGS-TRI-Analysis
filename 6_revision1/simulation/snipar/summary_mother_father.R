
setwd("~/work/family/simulation/snipar/08192025/")
path<-"~/work/family/simulation/snipar/08192025/"


files <- list.files(pattern = "^s.*\\.rda$")
true_mother = true_father = c(0,0,0,0)



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
  
  res_estimate_glm_fam_size1=array(0,c(sim_N,n_envir))
  res_sd_glm_fam_size1=array(0,c(sim_N,n_envir))
  res_p_glm_fam_size1=array(0,c(sim_N,n_envir))
  

  
  for(j in 1:sim_N){
    
    res_estimate_glm_fam_size1[j,]=res_glm_fam_parent[[j]]["G_mother",1]
    res_sd_glm_fam_size1[j,]=res_glm_fam_parent[[j]]["G_mother",2]
    res_p_glm_fam_size1[j,]=res_glm_fam_parent[[j]]["G_mother",4]
    
   
  }
  
  
  tmp=array(0,c(1,n_envir))
  
  tmp[1,]=apply(res_estimate_glm_fam_size1,2,mean)-true_mother[h]
  
  
  
  tmp=data.frame(tmp)
  colnames(tmp)[1]=c("bias_beta_prs")
  
  summary_bias[[h]]=tmp
  
  tmp=array(0,c(1,n_envir))
  
  tmp[1,]=apply(res_sd_glm_fam_size1,2,mean)
  
  tmp=data.frame(tmp)
  colnames(tmp)[1]=c("sd_beta_prs")
  
  summary_sd[[h]]=tmp
  
  tmp=array(0,c(1,n_envir))
  tmp[1,]=apply(res_p_glm_fam_size1,2,function(x) mean(x<0.05))
  
  
  
  tmp=data.frame(tmp)
  colnames(tmp)[1]=c("type1error_beta_prs")
  
  summary_type1error_power[[h]]=tmp
  
  tmp=array(0,c(1,n_envir))
  tmp[1,]=apply(res_estimate_glm_fam_size1,2,sd)
 
  
  tmp=data.frame(tmp)
  colnames(tmp)[1]=c("empirical_sd_beta_prs")
  
  summary_empirical_sd[[h]]=tmp
  
  tmp=array(0,c(1,n_envir))
  tmp[1,]=mean((res_estimate_glm_fam_size1-true_mother[h])^2)

  tmp=data.frame(tmp)
  colnames(tmp)[1]=c("mse_beta_prs")
  
  summary_mse[[h]]=tmp
  
  
  #95% coverage probability
  
  tmp=array(0,c(1,n_envir))
  upper=res_estimate_glm_fam_size1+qnorm(0.975)*res_sd_glm_fam_size1
  lower=res_estimate_glm_fam_size1-qnorm(0.975)*res_sd_glm_fam_size1
  tmp[1,1]=mean(ifelse(true_mother[h]<=upper[,1] & true_mother[h]>=lower[,1],1,0))

  tmp=data.frame(tmp)
  colnames(tmp)[1]=c("cp_beta_prs")
  
  summary_ci[[h]]=tmp
}


summary_mother=list()
tmp=data.frame(array(0,c(1,6)))
for(i in 1:4){
  tmp[,1]=summary_bias[[i]][,1]
  tmp[,2]=summary_sd[[i]][,1]
  tmp[,3]=summary_empirical_sd[[i]][,1]
  tmp[,4]=summary_type1error_power[[i]][,1]
  tmp[,5]=summary_ci[[i]][,1]
  tmp[,6]=summary_mse[[i]][,1]

  
  colnames(tmp)=c("bias","sd","empirical_sd","type1error/power","coverage_prob","mse")
  summary_mother[[i]]=data.frame(tmp)
  
}
summary_mother[[1]]$true_mother=0
summary_mother[[2]]$true_mother=0
summary_mother[[3]]$true_mother=0
summary_mother[[4]]$true_mother=0


summary_mother[[1]]$triosize=1000
summary_mother[[2]]$triosize=1000
summary_mother[[3]]$triosize=1000
summary_mother[[4]]$triosize=1000

summary_mother[[1]]$true_beta=0
summary_mother[[2]]$true_beta=0
summary_mother[[3]]$true_beta=0.4
summary_mother[[4]]$true_beta=0.4

summary_mother[[1]]$cor=0.4
summary_mother[[2]]$cor=0
summary_mother[[3]]$cor=0.4
summary_mother[[4]]$cor=0

summary_m = do.call(rbind,summary_mother)
save(summary_m,file=paste0(path,"results_summary_mother.rda"))






#Next check father
true_mother = true_father = c(0,0,0,0)



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
  
  res_estimate_glm_fam_size1=array(0,c(sim_N,n_envir))
  res_sd_glm_fam_size1=array(0,c(sim_N,n_envir))
  res_p_glm_fam_size1=array(0,c(sim_N,n_envir))
  
  
  
  for(j in 1:sim_N){
    
    res_estimate_glm_fam_size1[j,]=res_glm_fam_parent[[j]]["G_father",1]
    res_sd_glm_fam_size1[j,]=res_glm_fam_parent[[j]]["G_father",2]
    res_p_glm_fam_size1[j,]=res_glm_fam_parent[[j]]["G_father",4]
    
    
  }
  
  
  tmp=array(0,c(1,n_envir))
  
  tmp[1,]=apply(res_estimate_glm_fam_size1,2,mean)-true_father[h]
  
  
  
  tmp=data.frame(tmp)
  colnames(tmp)[1]=c("bias_beta_prs")
  
  summary_bias[[h]]=tmp
  
  tmp=array(0,c(1,n_envir))
  
  tmp[1,]=apply(res_sd_glm_fam_size1,2,mean)
  
  tmp=data.frame(tmp)
  colnames(tmp)[1]=c("sd_beta_prs")
  
  summary_sd[[h]]=tmp
  
  tmp=array(0,c(1,n_envir))
  tmp[1,]=apply(res_p_glm_fam_size1,2,function(x) mean(x<0.05))
  
  
  
  tmp=data.frame(tmp)
  colnames(tmp)[1]=c("type1error_beta_prs")
  
  summary_type1error_power[[h]]=tmp
  
  tmp=array(0,c(1,n_envir))
  tmp[1,]=apply(res_estimate_glm_fam_size1,2,sd)
  
  
  tmp=data.frame(tmp)
  colnames(tmp)[1]=c("empirical_sd_beta_prs")
  
  summary_empirical_sd[[h]]=tmp
  
  tmp=array(0,c(1,n_envir))
  tmp[1,]=mean((res_estimate_glm_fam_size1-true_father[h])^2)
  
  tmp=data.frame(tmp)
  colnames(tmp)[1]=c("mse_beta_prs")
  
  summary_mse[[h]]=tmp
  
  
  #95% coverage probability
  
  tmp=array(0,c(1,n_envir))
  upper=res_estimate_glm_fam_size1+qnorm(0.975)*res_sd_glm_fam_size1
  lower=res_estimate_glm_fam_size1-qnorm(0.975)*res_sd_glm_fam_size1
  tmp[1,1]=mean(ifelse(true_father[h]<=upper[,1] & true_father[h]>=lower[,1],1,0))
  
  tmp=data.frame(tmp)
  colnames(tmp)[1]=c("cp_beta_prs")
  
  summary_ci[[h]]=tmp
}


summary_father=list()
tmp=data.frame(array(0,c(1,6)))
for(i in 1:4){
  tmp[,1]=summary_bias[[i]][,1]
  tmp[,2]=summary_sd[[i]][,1]
  tmp[,3]=summary_empirical_sd[[i]][,1]
  tmp[,4]=summary_type1error_power[[i]][,1]
  tmp[,5]=summary_ci[[i]][,1]
  tmp[,6]=summary_mse[[i]][,1]
  
  
  colnames(tmp)=c("bias","sd","empirical_sd","type1error/power","coverage_prob","mse")
  summary_father[[i]]=data.frame(tmp)
  
}
summary_father[[1]]$true_father=0
summary_father[[2]]$true_father=0
summary_father[[3]]$true_father=0
summary_father[[4]]$true_father=0


summary_father[[1]]$triosize=1000
summary_father[[2]]$triosize=1000
summary_father[[3]]$triosize=1000
summary_father[[4]]$triosize=1000

summary_father[[1]]$true_beta=0
summary_father[[2]]$true_beta=0
summary_father[[3]]$true_beta=0.4
summary_father[[4]]$true_beta=0.4

summary_father[[1]]$cor=0.4
summary_father[[2]]$cor=0
summary_father[[3]]$cor=0.4
summary_father[[4]]$cor=0

summary_f = do.call(rbind,summary_father)
save(summary_f,file=paste0(path,"results_summary_father.rda"))


