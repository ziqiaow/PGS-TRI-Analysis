

#----------------------------------------------------------
#Figure S2
#May 21, 2024
#----------------------------------------------------------
setwd("./family/manuscript/plots/final/")
path<-"./family/manuscript/plots/final/"

library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyverse)
library(scales)
library(ggpubr)
library(cowplot)
library(reshape2)
library(RColorBrewer)

#####################################################################
#####################################################################
#####################################################################
# Panel A and B

#Power for PGS main effect

filenames <- list.files("./family/simulation/power/redo_06112024/", pattern=".*power.*\\.rda$", full.names=TRUE)
sim_N=1000
power_beta=array(0,c(length(filenames),4))
power_beta_glm=array(0,c(length(filenames),4))
power_beta_tdt=array(0,c(length(filenames),4))
res_p_fam=array(0,c(sim_N,4))
res_p_tdt=array(0,c(sim_N,4))
res_p_glm=array(0,c(sim_N,4))
for(i in 1:length(filenames)){
  
  load(filenames[i])
  for(j in 1:sim_N){
    
    res_p_fam[j,1]=res_fam_size1[[j]][,4]
    res_p_fam[j,2]=res_fam_size2[[j]][,4]
    res_p_fam[j,3]=res_fam_size3[[j]][,4]
    res_p_fam[j,4]=res_fam_size4[[j]][,4]
    
    
    res_p_glm[j,1]=res_glm_fam_size1[[j]][2,4]
    res_p_glm[j,2]=res_glm_fam_size2[[j]][2,4]
    res_p_glm[j,3]=res_glm_fam_size3[[j]][2,4]
    res_p_glm[j,4]=res_glm_fam_size4[[j]][2,4]
    
    res_p_tdt[j,1]=res_ptdt_size1[[j]][,4]
    res_p_tdt[j,2]=res_ptdt_size2[[j]][,4]
    res_p_tdt[j,3]=res_ptdt_size3[[j]][,4]
    res_p_tdt[j,4]=res_ptdt_size4[[j]][,4]
    
  }
  
  power_beta[i,]=apply(res_p_fam,2,function(x) mean(x<0.05))
  power_beta_glm[i,]=apply(res_p_glm,2,function(x) mean(x<0.05))
  power_beta_tdt[i,]=apply(res_p_tdt,2,function(x) mean(x<0.05))
}



power_plot=data.frame(power_beta)
colnames(power_plot)=c("200","500","1000","2000")
power_plot$logor=seq(0,0.7,0.1)

power_plot_tdt=data.frame(power_beta_tdt)
colnames(power_plot_tdt)=c("200","500","1000","2000")
power_plot_tdt$logor=seq(0,0.7,0.1)

power_plot_glm=data.frame(power_beta_glm)
colnames(power_plot_glm)=c("200","500","1000","2000")
power_plot_glm$logor=seq(0,0.7,0.1)

power_plot=rbind(power_plot,power_plot_tdt,power_plot_glm)
power_plot$Method=c(rep("PGS-TRI",8),rep("pTDT",8),rep("Logistic Regression",8))
power_plot$Method=factor(power_plot$Method,level=c("PGS-TRI","Logistic Regression","pTDT"))

# Specify id.vars: the variables to keep but not split apart on
power_plot<-melt(power_plot, id.vars=c("logor","Method"))
colnames(power_plot)[3]=c("trio")
p1 <- ggplot(power_plot[which(power_plot$trio=="1000" & power_plot$logor<=0.5) ,],aes(logor,value,color = Method))+ylim(c(0,1.1))  +geom_line(size=2,position = position_dodge(width = 0.005))+
  geom_point(size=3) + scale_x_continuous('True log OR',breaks=seq(0,0.7,0.1)) +
  scale_y_continuous('Power',breaks=c(0.05,0.5,0.8,1))+ labs(color='Method')  +theme_bw()+ ggtitle(expression(bold(paste("\u03B2(PGS)"))))+                                                               # Change font size
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "orange", size=1.2)+
scale_color_manual(values = c(pal_nejm("default")(5)[1],"Royal Blue",pal_nejm("default")(5)[c(4)]))+
theme(plot.title = element_text(color = "black", size = 14, face = "bold")) +  theme(strip.text.x = element_text(face = "bold", size = 14),strip.text.y = element_text(face = "bold", size = 14))+theme(legend.position = c(0.7, 0.25))+  theme(legend.text = element_text(face = "bold", size=14),legend.title = element_text(face = "bold", size=14))+theme(axis.title.x = element_text(face="bold",size=14, color="black"), axis.text.x = element_text(face="bold",size=14, color="black"), axis.title.y = element_text(face="bold",size=14, color="black"), axis.text.y = element_text(face="bold", size=14, color="black"))


# Data plotting
library(reshape2)
library(ggplot2)
library("ggsci")
power_plot=data.frame(power_beta)
colnames(power_plot)=c("200","500","1000","2000")
power_plot$logor=seq(0,0.7,0.1)

# Specify id.vars: the variables to keep but not split apart on
power_plot<-melt(power_plot, id.vars=c("logor"))
colnames(power_plot)[2]=c("trio")
cols <- brewer.pal(9,"YlOrRd")
p2 <- ggplot(power_plot,aes(logor,value,color = trio))+ylim(c(0,1.1)) +geom_line(size=1.2) +
  geom_point(size=3) + scale_x_continuous('True log OR',breaks=seq(0,0.7,0.1)) +
  scale_y_continuous('Power',breaks=c(0.05,0.5,0.8,1))+labs(color='# of Trios')  + ggtitle(expression(paste("\u03B2(PGS)")))+   # geom_hline(yintercept=0.8, linetype="dashed",  color = "orange", size=1.2)
  theme_bw()+                                                            # Change font size
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "orange", size=1.2)+ #+ scale_color_npg()
  scale_color_manual(values = cols[4:7])+theme(strip.text.x = element_text(face = "bold", size = 12),strip.text.y = element_text(face = "bold", size = 12))+theme(legend.position=c(0.9,0.2))+  theme(legend.text = element_text(face = "bold", size=12),legend.title = element_text(face = "bold", size=12))+theme(axis.title.x = element_text(face="bold",size=12, color="black"), axis.text.x = element_text(face="bold",size=12, color="black"), axis.title.y = element_text(face="bold",size=12, color="black"), axis.text.y = element_text(face="bold", size=10, color="black"))


#####################################################################
#####################################################################
#####################################################################
# Panel C and D


##Power for Differential Indirect Genetic effect

sim_N=1000
power_delta_glm=array(0,c(8,4))
power_delta=array(0,c(8,4))
res_p_glm=array(0,c(sim_N,4))
res_p_fam_delta=array(0,c(sim_N,4))

for(i in 1:8){
  load(paste0("./family/simulation/nurture/redo/05212024/power/scenario",i,".rda"))
  
  for(j in 1:sim_N){
    
    
    res_p_fam_delta[j,1]=res_fam_size1[[j]]$res_delta[,4]
    res_p_fam_delta[j,2]=res_fam_size2[[j]]$res_delta[,4]
    res_p_fam_delta[j,3]=res_fam_size3[[j]]$res_delta[,4]
    res_p_fam_delta[j,4]=res_fam_size4[[j]]$res_delta[,4]
    
    res_p_glm[j,1]=res_glm_fam_size1[[j]]$res_delta[,4]
    res_p_glm[j,2]=res_glm_fam_size2[[j]]$res_delta[,4]
    res_p_glm[j,3]=res_glm_fam_size3[[j]]$res_delta[,4]
    res_p_glm[j,4]=res_glm_fam_size4[[j]]$res_delta[,4]
    
  }
  
  
  power_delta[i,]=apply(res_p_fam_delta,2,function(x) mean(x<0.05))
  power_delta_glm[i,]=apply(res_p_glm,2,function(x) mean(x<0.05))
  
  
}


power_plot=data.frame(power_delta)
colnames(power_plot)=c("200","500","1000","2000")
power_plot$logor=seq(0,0.7,0.1)


power_plot_glm=data.frame(power_delta_glm)
colnames(power_plot_glm)=c("200","500","1000","2000")
power_plot_glm$logor=seq(0,0.7,0.1)

power_plot=rbind(power_plot,power_plot_glm)
power_plot$Method=c(rep("PGS-TRI",8),rep("Logistic Regression",8))
power_plot$Method=factor(power_plot$Method,level=c("PGS-TRI","Logistic Regression"))

# Specify id.vars: the variables to keep but not split apart on
power_plot<-melt(power_plot, id.vars=c("logor","Method"))
colnames(power_plot)[3]=c("trio")
p3 <- ggplot(power_plot[which(power_plot$trio=="1000" & power_plot$logor<=0.5) ,],aes(logor,value,color = Method))+ylim(c(0,1.1))  +geom_line(size=2,position = position_dodge(width = 0.005))+
  geom_point(size=3) + scale_x_continuous('True log OR',breaks=seq(0,0.7,0.1)) +
  scale_y_continuous('Power',breaks=c(0.05,0.5,0.8,1))+ labs(color='Method') +theme_bw()+ ggtitle(expression(bold(paste("\u03B4-IDE"))))+                                                               # Change font size
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "orange", size=1.2)+theme(plot.title = element_text(color = "black", size = 14, face = "bold")) + 
  scale_color_manual(values = c(pal_nejm("default")(5)[1],"Royal Blue",pal_nejm("default")(5)[c(4)]))+  theme(strip.text.x = element_text(face = "bold", size = 14),strip.text.y = element_text(face = "bold", size = 14))+theme(legend.position = c(0.7, 0.25))+  theme(legend.text = element_text(face = "bold", size=14),legend.title = element_text(face = "bold", size=14))+theme(axis.title.x = element_text(face="bold",size=14, color="black"), axis.text.x = element_text(face="bold",size=14, color="black"), axis.title.y = element_text(face="bold",size=14, color="black"), axis.text.y = element_text(face="bold", size=14, color="black"))



power_plot=data.frame(power_delta)
colnames(power_plot)=c("200","500","1000","2000")
power_plot$logor=seq(0,0.7,0.1)

# Specify id.vars: the variables to keep but not split apart on
power_plot<-melt(power_plot, id.vars=c("logor"))
colnames(power_plot)[2]=c("trio")
p4 <- ggplot(power_plot,aes(logor,value,color = trio))+ylim(c(0,1.1)) +geom_line(size=1.2) +
  geom_point(size=3) + scale_x_continuous('True log OR',breaks=seq(0,0.7,0.1)) +
  scale_y_continuous('Power',breaks=c(0.05,0.5,0.8,1))+labs(color='# of Trios')  + ggtitle(expression(paste("\u03B4-IDE")))+   # geom_hline(yintercept=0.8, linetype="dashed",  color = "orange", size=1.2)
  theme_bw()+                                                            # Change font size
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "orange", size=1.2)+ #+ scale_color_npg()
  scale_color_manual(values = cols[4:7])+theme(strip.text.x = element_text(face = "bold", size = 12),strip.text.y = element_text(face = "bold", size = 12))+theme(legend.position=c(0.9,0.2))+  theme(legend.text = element_text(face = "bold", size=12),legend.title = element_text(face = "bold", size=12))+theme(axis.title.x = element_text(face="bold",size=12, color="black"), axis.text.x = element_text(face="bold",size=12, color="black"), axis.title.y = element_text(face="bold",size=12, color="black"), axis.text.y = element_text(face="bold", size=10, color="black"))


#####################################################################
#####################################################################
#####################################################################
# Panel E and F

#Power for PGS x E1



sim_N=1000
power_gxe_glm=array(0,c(8,4))
power_gxe=array(0,c(8,4))
power_caseonly=array(0,c(8,4))
res_p_glm=array(0,c(sim_N,4))
res_p_fam_gxe=array(0,c(sim_N,4))
res_p_caseonly=array(0,c(sim_N,4))

for(i in 1:8){
  load(paste0("./family/simulation_gxe/07232024/power/e1/scenario",i,".rda"))
  
  for(j in 1:sim_N){
    
    
    res_p_fam_gxe[j,1]=res_fam_size1[[j]]$res_beta[2,4]
    res_p_fam_gxe[j,2]=res_fam_size2[[j]]$res_beta[2,4]
    res_p_fam_gxe[j,3]=res_fam_size3[[j]]$res_beta[2,4]
    res_p_fam_gxe[j,4]=res_fam_size4[[j]]$res_beta[2,4]
    
    res_p_glm[j,1]=res_glm_fam_size1[[j]][4,4]
    res_p_glm[j,2]=res_glm_fam_size2[[j]][4,4]
    res_p_glm[j,3]=res_glm_fam_size3[[j]][4,4]
    res_p_glm[j,4]=res_glm_fam_size4[[j]][4,4]
    
    
    res_p_caseonly[j,1]=res_caseonly_size1[[j]][2,4]
    res_p_caseonly[j,2]=res_caseonly_size2[[j]][2,4]
    res_p_caseonly[j,3]=res_caseonly_size3[[j]][2,4]
    res_p_caseonly[j,4]=res_caseonly_size4[[j]][2,4]
    
  }
  
  
  power_gxe[i,]=apply(res_p_fam_gxe,2,function(x) mean(x<0.05))
  power_gxe_glm[i,]=apply(res_p_glm,2,function(x) mean(x<0.05))
  power_caseonly[i,]=apply(res_p_caseonly,2,function(x) mean(x<0.05))
  
  
  
}




# Data plotting
library(reshape2)
library(ggplot2)
library("ggsci")
power_plot=data.frame(power_gxe)
colnames(power_plot)=c("200","500","1000","2000")
power_plot$logor=seq(0,0.7,0.1)


power_plot_glm=data.frame(power_gxe_glm)
colnames(power_plot_glm)=c("200","500","1000","2000")
power_plot_glm$logor=seq(0,0.7,0.1)


power_plot_caseonly=data.frame(power_caseonly)
colnames(power_plot_caseonly)=c("200","500","1000","2000")
power_plot_caseonly$logor=seq(0,0.7,0.1)


power_plot=rbind(power_plot,power_plot_glm,power_plot_caseonly)
power_plot$Method=c(rep("PGS-TRI",8),rep("Logistic Regression",8),rep("Case-Only",8))
power_plot$Method=factor(power_plot$Method,level=c("PGS-TRI","Logistic Regression","Case-Only"))

# Specify id.vars: the variables to keep but not split apart on
power_plot<-melt(power_plot, id.vars=c("logor","Method"))
colnames(power_plot)[3]=c("trio")

p5 <- ggplot(power_plot[which(power_plot$trio=="1000" & power_plot$logor<=0.5) ,],aes(logor,value,color = Method))+ylim(c(0,1.1))  +geom_line(size=2,position = position_dodge(width = 0.005))+
  geom_point(size=3) + scale_x_continuous('True log OR',breaks=seq(0,0.7,0.1)) +
  scale_y_continuous('Power',breaks=c(0.05,0.5,0.8,1))+ labs(color='Method') +theme_bw()+ ggtitle(expression(bold( paste("\u03B2(PGSxE1)"))))+                                                               # Change font size
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "orange", size=1.2)+theme(plot.title = element_text(color = "black", size = 14, face = "bold")) +      
  scale_color_manual(values = c(pal_nejm("default")(5)[1],"Royal Blue",pal_jco("default")(3)[2]))+theme(strip.text.x = element_text(face = "bold", size = 14),strip.text.y = element_text(face = "bold", size = 14))+theme(legend.position = c(0.7, 0.25))+  theme(legend.text = element_text(face = "bold", size=14),legend.title = element_text(face = "bold", size=14))+theme(axis.title.x = element_text(face="bold",size=14, color="black"), axis.text.x = element_text(face="bold",size=14, color="black"), axis.title.y = element_text(face="bold",size=14, color="black"), axis.text.y = element_text(face="bold", size=14, color="black"))


# Data plotting
library(reshape2)
library(ggplot2)
library("ggsci")
power_plot=data.frame(power_gxe)
colnames(power_plot)=c("200","500","1000","2000")
power_plot$logor=seq(0,0.7,0.1)

# Specify id.vars: the variables to keep but not split apart on
power_plot<-melt(power_plot, id.vars=c("logor"))
colnames(power_plot)[2]=c("trio")

p6 <- ggplot(power_plot,aes(logor,value,color = trio))+ylim(c(0,1.1)) +geom_line(size=1.2) +
  geom_point(size=3) + scale_x_continuous('True log OR',breaks=seq(0,0.7,0.1)) +
  scale_y_continuous('Power',breaks=c(0.05,0.5,0.8,1))+labs(color='# of Trios') + ggtitle(expression(paste("\u03B2(PGSxE1)")))+   # geom_hline(yintercept=0.8, linetype="dashed",  color = "orange", size=1.2)
  theme_bw()+                                                            # Change font size
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "orange", size=1.2)+ #+ scale_color_npg()
  scale_color_manual(values = cols[4:7])+theme(strip.text.x = element_text(face = "bold", size = 12),strip.text.y = element_text(face = "bold", size = 12))+theme(legend.position=c(0.9,0.2))+  theme(legend.text = element_text(face = "bold", size=12),legend.title = element_text(face = "bold", size=12))+theme(axis.title.x = element_text(face="bold",size=12, color="black"), axis.text.x = element_text(face="bold",size=12, color="black"), axis.title.y = element_text(face="bold",size=12, color="black"), axis.text.y = element_text(face="bold", size=10, color="black"))




#####################################################################
#####################################################################
#####################################################################
# Panel G and H

#Power for PGS x E2



sim_N=1000
power_gxe_glm=array(0,c(8,4))
power_gxe=array(0,c(8,4))
power_caseonly=array(0,c(8,4))
res_p_glm=array(0,c(sim_N,4))
res_p_fam_gxe=array(0,c(sim_N,4))
res_p_caseonly=array(0,c(sim_N,4))

for(i in 1:8){
  load(paste0("./family/simulation_gxe/07232024/power/e2/scenario",i,".rda"))
  
  for(j in 1:sim_N){
    
    
    res_p_fam_gxe[j,1]=res_fam_size1[[j]]$res_beta[2,4]
    res_p_fam_gxe[j,2]=res_fam_size2[[j]]$res_beta[2,4]
    res_p_fam_gxe[j,3]=res_fam_size3[[j]]$res_beta[2,4]
    res_p_fam_gxe[j,4]=res_fam_size4[[j]]$res_beta[2,4]
    
    res_p_glm[j,1]=res_glm_fam_size1[[j]][4,4]
    res_p_glm[j,2]=res_glm_fam_size2[[j]][4,4]
    res_p_glm[j,3]=res_glm_fam_size3[[j]][4,4]
    res_p_glm[j,4]=res_glm_fam_size4[[j]][4,4]
    
    
    res_p_caseonly[j,1]=res_caseonly_size1[[j]][2,4]
    res_p_caseonly[j,2]=res_caseonly_size2[[j]][2,4]
    res_p_caseonly[j,3]=res_caseonly_size3[[j]][2,4]
    res_p_caseonly[j,4]=res_caseonly_size4[[j]][2,4]
    
  }
  
  
  power_gxe[i,]=apply(res_p_fam_gxe,2,function(x) mean(x<0.05))
  power_gxe_glm[i,]=apply(res_p_glm,2,function(x) mean(x<0.05))
  power_caseonly[i,]=apply(res_p_caseonly,2,function(x) mean(x<0.05))
  
  
  
}




power_plot=data.frame(power_gxe)
colnames(power_plot)=c("200","500","1000","2000")
#rownames(power_plot)=seq(0.1,0.7,0.1)
power_plot$logor=seq(0,0.7,0.1)


power_plot_glm=data.frame(power_gxe_glm)
colnames(power_plot_glm)=c("200","500","1000","2000")
#rownames(power_plot)=seq(0.1,0.7,0.1)
power_plot_glm$logor=seq(0,0.7,0.1)


power_plot_caseonly=data.frame(power_caseonly)
colnames(power_plot_caseonly)=c("200","500","1000","2000")
#rownames(power_plot)=seq(0.1,0.7,0.1)
power_plot_caseonly$logor=seq(0,0.7,0.1)


power_plot=rbind(power_plot,power_plot_glm,power_plot_caseonly)
power_plot$Method=c(rep("PGS-TRI",8),rep("Logistic Regression",8),rep("Case-Only",8))
power_plot$Method=factor(power_plot$Method,level=c("PGS-TRI","Logistic Regression","Case-Only"))

# Specify id.vars: the variables to keep but not split apart on
power_plot<-melt(power_plot, id.vars=c("logor","Method"))
colnames(power_plot)[3]=c("trio")

p7 <- ggplot(power_plot[which(power_plot$trio=="1000" & power_plot$logor<=0.5) ,],aes(logor,value,color = Method))+ylim(c(0,1.1))  +geom_line(size=2,position = position_dodge(width = 0.005))+
  geom_point(size=3) + scale_x_continuous('True log OR',breaks=seq(0,0.7,0.1)) +
  scale_y_continuous('Power',breaks=c(0.05,0.5,0.8,1))+ labs(color='Method')  +theme_bw()+ ggtitle(expression(bold( paste("\u03B2(PGSxE2)")) ))+                                                               # Change font size
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "orange", size=1.2)+ theme(plot.title = element_text(color = "black", size = 14, face = "bold")) +      
  scale_color_manual(values = c(pal_nejm("default")(5)[1],"Royal Blue",pal_jco("default")(3)[2]))+theme(strip.text.x = element_text(face = "bold", size = 14),strip.text.y = element_text(face = "bold", size = 14))+theme(legend.position = c(0.7, 0.25))+  theme(legend.text = element_text(face = "bold", size=14),legend.title = element_text(face = "bold", size=14))+theme(axis.title.x = element_text(face="bold",size=14, color="black"), axis.text.x = element_text(face="bold",size=14, color="black"), axis.title.y = element_text(face="bold",size=14, color="black"), axis.text.y = element_text(face="bold", size=14, color="black"))



power_plot=data.frame(power_gxe)
colnames(power_plot)=c("200","500","1000","2000")
#rownames(power_plot)=seq(0.1,0.7,0.1)
power_plot$logor=seq(0,0.7,0.1)

# Specify id.vars: the variables to keep but not split apart on
power_plot<-melt(power_plot, id.vars=c("logor"))
colnames(power_plot)[2]=c("trio")

p8 <- ggplot(power_plot,aes(logor,value,color = trio))+ylim(c(0,1.1)) +geom_line(size=1.2) +
  geom_point(size=3) + scale_x_continuous('True log OR',breaks=seq(0,0.7,0.1)) +
  scale_y_continuous('Power',breaks=c(0.05,0.5,0.8,1))+labs(color='# of Trios') + ggtitle(expression(paste("\u03B2(PGSxE2)")))+   # geom_hline(yintercept=0.8, linetype="dashed",  color = "orange", size=1.2)
  theme_bw()+                                                            # Change font size
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "orange", size=1.2)+ #+ scale_color_npg()
  scale_color_manual(values = cols[4:7])+theme(strip.text.x = element_text(face = "bold", size = 12),strip.text.y = element_text(face = "bold", size = 12))+theme(legend.position=c(0.9,0.2))+  theme(legend.text = element_text(face = "bold", size=12),legend.title = element_text(face = "bold", size=12))+theme(axis.title.x = element_text(face="bold",size=12, color="black"), axis.text.x = element_text(face="bold",size=12, color="black"), axis.title.y = element_text(face="bold",size=12, color="black"), axis.text.y = element_text(face="bold", size=10, color="black"))





#####################################################################
#####################################################################
#####################################################################
# legends
power_plot$Method = rep(c("PGS-TRI","Logistic Regression","pTDT","Case-Only"),8)
power_plot$Method=factor(power_plot$Method,level=c("PGS-TRI","Logistic Regression","pTDT","Case-Only"))


p_tmp <- ggplot(power_plot[which(power_plot$trio=="1000" & power_plot$logor<=0.5) ,],aes(logor,value,color = Method))+ylim(c(0,1.1))  +geom_line(size=2,position = position_dodge(width = 0.005))+
  geom_point(size=3) + scale_x_continuous('True log OR',breaks=seq(0,0.7,0.1)) +
  scale_y_continuous('Power',breaks=c(0.05,0.5,0.8,1))+ labs(color='Method')  +theme_bw()+ ggtitle(expression(paste("\u03B2(PGSxE2)")))+                                                               # Change font size
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "orange", size=1.2)+
  scale_color_manual(values = c(pal_nejm("default")(5)[1],"Royal Blue",pal_nejm("default")(5)[c(4)],pal_jco("default")(3)[2]))+theme(strip.text.x = element_text(face = "bold", size = 12),strip.text.y = element_text(face = "bold", size = 12))+theme(legend.position = c(0.7, 0.25))+  theme(legend.text = element_text(face = "bold", size=12),legend.title = element_text(face = "bold", size=12))+theme(axis.title.x = element_text(face="bold",size=12, color="black"), axis.text.x = element_text(face="bold",size=12, color="black"), axis.title.y = element_text(face="bold",size=12, color="black"), axis.text.y = element_text(face="bold", size=10, color="black"))


p_tmp <- p_tmp+  theme(strip.text.x = element_text(face = "bold", size = 18),strip.text.y = element_text(face = "bold", size = 18))+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size = 18),legend.title = element_text(face = "bold", size = 18))+theme(axis.title.x = element_text(face="bold",size = 18, color="black"), axis.text.x = element_text(face="bold",size = 18, color="black",angle = 45, vjust = 0.5, hjust = 0.5), axis.title.y = element_text(face="bold",size = 18, color="black"), axis.text.y = element_text(face="bold", size = 18, color="black"))

legend_b <- cowplot::get_plot_component(p_tmp+ theme(legend.position="top"), 'guide-box-top', return_all = TRUE)

#####################################################################
#####################################################################
#####################################################################


#Finally, put the plots together


p1 <- p1+ theme(legend.position="none")
#p2 <- p2+ theme(legend.position="none")
p3 <- p3+ theme(legend.position="none")
#p4 <- p4+ theme(legend.position="none")
p5 <- p5+ theme(legend.position="none")
#p6 <- p6+ theme(legend.position="none")
p7 <- p7+ theme(legend.position="none")



p <- ggarrange(ggarrange(p1,p3,
                         ncol = 2, labels = c("a", "b"),
                         widths = c(0.5,0.5)) , 
               ggarrange(p5,p7,
                         ncol = 2, labels = c("c","d"),
                         widths = c(0.5,0.5)) ,
               nrow = 2, 
               labels = c(NA, NA),
               heights = c(0.5,0.5))

p <- plot_grid( legend_b,p,  ncol = 1, rel_heights = c(.05, 1))

ggsave(filename=paste0("Extended_Figure2_09262024.pdf"),
       plot=p, 
       path="./family/manuscript/plots/final/",
       width=380, height=300, units="mm", dpi=320,device = cairo_pdf)

ggsave(filename=paste0("Extended_Figure2_09262024.png"),
       plot=p, 
       path="./family/manuscript/plots/final/",
       width=380, height=300, units="mm", dpi=320)

