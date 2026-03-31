

#----------------------------------------------------------
#Figure S2
#----------------------------------------------------------
setwd("./family/manuscript/draft4/plots/extended/")
path<-"./family/manuscript/draft4/plots/extended/"

library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyverse)
library(scales)
library(ggpubr)
library(cowplot)
library(reshape2)
library(RColorBrewer)



My_Theme = theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  title = element_text(size = 7),
  text = element_text(color = "black",size = 7),
  axis.text = element_text(color = "black", size = 7,face = "bold"),
  axis.title = element_text(size = 7,face = "bold"),
  # axis.text.x = element_text(angle = 45,  vjust = 0.55, hjust = 0.5),
  strip.text.x = element_text(size = 7,face = "bold"),strip.text.y = element_text(size = 7,face = "bold"),
  panel.border = element_rect(size = 0.15),
  strip.background = element_rect(size = 0.15),
  axis.line = element_line(size = 0.15),
  axis.ticks = element_line(size = 0.15)
  
)

#####################################################################
#####################################################################
#####################################################################
# Panel A

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
#rownames(power_plot)=seq(0.1,0.7,0.1)
power_plot$logor=seq(0,0.7,0.1)

power_plot_tdt=data.frame(power_beta_tdt)
colnames(power_plot_tdt)=c("200","500","1000","2000")
#rownames(power_plot)=seq(0.1,0.7,0.1)
power_plot_tdt$logor=seq(0,0.7,0.1)

power_plot_glm=data.frame(power_beta_glm)
colnames(power_plot_glm)=c("200","500","1000","2000")
#rownames(power_plot)=seq(0.1,0.7,0.1)
power_plot_glm$logor=seq(0,0.7,0.1)

power_plot=rbind(power_plot,power_plot_tdt,power_plot_glm)
power_plot$Method=c(rep("PGS-TRI",8),rep("pTDT",8),rep("Logistic Regression",8))
power_plot$Method=factor(power_plot$Method,level=c("PGS-TRI","Logistic Regression","pTDT"))

# Specify id.vars: the variables to keep but not split apart on
power_plot<-melt(power_plot, id.vars=c("logor","Method"))
colnames(power_plot)[3]=c("trio")
p1 <- ggplot(power_plot[which(power_plot$trio=="1000" & power_plot$logor<=0.5) ,],aes(logor,value,color = Method))+ylim(c(0,1.1))  +geom_line(size=1,position = position_dodge(width = 0.005))+
  geom_point(size=1.5) + scale_x_continuous('True log OR',breaks=seq(0,0.7,0.1)) +
  scale_y_continuous('Power',breaks=c(0.05,0.5,0.8,1))+ labs(color='Method')  +theme_bw()+ ggtitle(expression(bold(paste("\u03B2(PGS)"))))+                                                               # Change font size
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "orange", size=1)+
scale_color_manual(values = c(pal_nejm("default")(5)[1],"Royal Blue",pal_nejm("default")(5)[c(4)]))+ theme_bw()+      
  theme(plot.title = element_text(color = "black", size = 7, face = "bold")) +                                                          
  My_Theme+theme(legend.position="") 



#####################################################################
#####################################################################
#####################################################################
# Panel B


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
p3 <- ggplot(power_plot[which(power_plot$trio=="1000" & power_plot$logor<=0.5) ,],aes(logor,value,color = Method))+ylim(c(0,1.1))  +geom_line(size=1,position = position_dodge(width = 0.005))+
  geom_point(size=1.5) + scale_x_continuous('True log OR',breaks=seq(0,0.7,0.1)) +
  scale_y_continuous('Power',breaks=c(0.05,0.5,0.8,1))+ labs(color='Method') +theme_bw()+ ggtitle(expression(bold(paste("\u03B4-IDE"))))+                                                               # Change font size
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "orange", size=1)+
  scale_color_manual(values = c(pal_nejm("default")(5)[1],"Royal Blue",pal_nejm("default")(5)[c(4)]))+ theme_bw()+      
  theme(plot.title = element_text(color = "black", size = 7, face = "bold")) +                                                          
  My_Theme+theme(legend.position="") 




#####################################################################
#####################################################################
#####################################################################
# Panel C

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

p5 <- ggplot(power_plot[which(power_plot$trio=="1000" & power_plot$logor<=0.5) ,],aes(logor,value,color = Method))+ylim(c(0,1.1))  +geom_line(size=1,position = position_dodge(width = 0.005))+
  geom_point(size=1.5) + scale_x_continuous('True log OR',breaks=seq(0,0.7,0.1)) +
  scale_y_continuous('Power',breaks=c(0.05,0.5,0.8,1))+ labs(color='Method') +theme_bw()+ ggtitle(expression(bold( paste("\u03B2(PGSxE1)"))))+                                                               # Change font size
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "orange", size=1) +      
  scale_color_manual(values = c(pal_nejm("default")(5)[1],"Royal Blue",pal_jco("default")(3)[2]))+ theme_bw()+      
  theme(plot.title = element_text(color = "black", size = 7, face = "bold")) +                                                          
  My_Theme+theme(legend.position="") 




#####################################################################
#####################################################################
#####################################################################
# Panel D

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

p7 <- ggplot(power_plot[which(power_plot$trio=="1000" & power_plot$logor<=0.5) ,],aes(logor,value,color = Method))+ylim(c(0,1.1))  +geom_line(size=1,position = position_dodge(width = 0.005))+
  geom_point(size=1.5) + scale_x_continuous('True log OR',breaks=seq(0,0.7,0.1)) +
  scale_y_continuous('Power',breaks=c(0.05,0.5,0.8,1))+ labs(color='Method')  +theme_bw()+ ggtitle(expression(bold( paste("\u03B2(PGSxE2)")) ))+                                                               # Change font size
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "orange", size=1)+     
  scale_color_manual(values = c(pal_nejm("default")(5)[1],"Royal Blue",pal_jco("default")(3)[2]))+ theme_bw()+      
  theme(plot.title = element_text(color = "black", size = 7, face = "bold")) +                                                          
  My_Theme+theme(legend.position="") 



power_plot=data.frame(power_gxe)
colnames(power_plot)=c("200","500","1000","2000")
#rownames(power_plot)=seq(0.1,0.7,0.1)
power_plot$logor=seq(0,0.7,0.1)

# Specify id.vars: the variables to keep but not split apart on
power_plot<-melt(power_plot, id.vars=c("logor"))
colnames(power_plot)[2]=c("trio")





#####################################################################
#####################################################################
#####################################################################
# legends
power_plot$Method = rep(c("PGS-TRI","Logistic Regression","pTDT","Case-Only"),8)
power_plot$Method=factor(power_plot$Method,level=c("PGS-TRI","Logistic Regression","pTDT","Case-Only"))


p_tmp <- ggplot(power_plot[which(power_plot$trio=="1000" & power_plot$logor<=0.5) ,],aes(logor,value,color = Method))+ylim(c(0,1.1))  +geom_line(size=1,position = position_dodge(width = 0.005))+
  geom_point(size=1.5) + scale_x_continuous('True log OR',breaks=seq(0,0.7,0.1)) +
  scale_y_continuous('Power',breaks=c(0.05,0.5,0.8,1))+ labs(color='Method')  +theme_bw()+ ggtitle(expression(paste("\u03B2(PGSxE2)")))+                                                               # Change font size
  geom_hline(yintercept=0.05, linetype="dashed", 
             color = "orange", size=1)+
  scale_color_manual(values = c(pal_nejm("default")(5)[1],"Royal Blue",pal_nejm("default")(5)[c(4)],pal_jco("default")(3)[2]))

p_tmp <- p_tmp+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size = 7),legend.title = element_text(face = "bold", size = 7),legend.key.size = unit(3, 'mm'))

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
                         ncol = 2, labels = c("a", "b"),font.label = list(size = 7),
                         widths = c(0.5,0.5)) , 
               ggarrange(p5,p7,
                         ncol = 2, labels = c("c","d"),font.label = list(size = 7),
                         widths = c(0.5,0.5)) ,
               nrow = 2, 
               labels = c(NA, NA),font.label = list(size = 7),
               heights = c(0.5,0.5))

p <- plot_grid( legend_b,p,  ncol = 1, rel_heights = c(.05, 1))

ggsave(filename=paste0("Extended_Figure2.pdf"),
       plot=p, 
       path="./family/manuscript/draft4/plots/extended/",
       width=180, height=160, units="mm", dpi=400,device = cairo_pdf)

ggsave(filename=paste0("Extended_Figure2.png"),
       plot=p, 
       path="./family/manuscript/draft4/plots/extended/",
       width=180, height=160, units="mm", dpi=400)

