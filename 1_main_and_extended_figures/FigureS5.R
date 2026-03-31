#----------------------------------------------------------
#Extended Figure 5
#----------------------------------------------------------
setwd("./family/manuscript/plots/draft2/")
path<-"./family/manuscript/plots/draft2/"

library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyverse)
library(scales)
library(ggpubr)
library(cowplot)


My_Theme = theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  title = element_text(size = 7),
  text = element_text(color = "black",size = 7),
  axis.text = element_text(color = "black", size = 7,face = "bold"),
  axis.title = element_text(size = 7,face = "bold"),
  axis.text.x = element_text(angle = 45,  vjust = 0.55, hjust = 0.5),
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

#type I error for UKB simulation results for DE
#load data with PS matched by geographical regions only

load("./family/simulation/snipar/08192025/results_summary_beta_main.rda")

#type I error
data_plot= bind_rows(summary_beta[1])#bind_rows(summary_beta[1:2])
data_plot = data_plot[which(data_plot$method %in% c("logistic","ptdt","fam","logistic_parent_adj","fam_parent_adj")),]

data_plot$method = c("PGS-TRI (DE only)","Logistic Regression","PGS-TRI (DE+IDE)","Logistic Regression + PGS(MF)","pTDT")
data_plot$method=factor(data_plot$method,levels=c("PGS-TRI (DE only)","PGS-TRI (DE+IDE)","Logistic Regression","Logistic Regression + PGS(MF)","pTDT"))


colnames(data_plot)[7]="Method"

data_plot$effect = "DE"
data_plot$label = "\u03B2(G)=0,\u03B4-IDE=0"

p1 <- ggplot(data=data_plot, aes(x=Method, y=type1error.power , fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width = 0.9,linewidth = 0.1)+
  scale_fill_manual(values = c(pal_nejm("default")(5)[c(1,3)],"Royal Blue","Light Sky Blue",pal_nejm("default")(5)[c(4)])) +
  scale_x_discrete(NULL) +
  scale_y_continuous('Type I Error for DE',breaks = c(0.05,0.2),limits = c(0,0.2))+
  geom_hline(yintercept=0.05, col = "coral",linetype = "dashed",size=0.5)  + 
  facet_grid(effect~label) + theme_bw()+      
  theme(plot.title = element_text(color = "black", size = 7, face = "bold")) +                                                          
  My_Theme+theme(legend.position="") 
 



# Panel B

#type I error for UKB simulation results for IDE
#load data with PS matched by geographical regions only

load("./family/simulation/snipar/08192025/results_summary_delta_indirect.rda")
load("./family/simulation/snipar/08192025/results_summary_mother.rda")
load("./family/simulation/snipar/08192025/results_summary_father.rda")

#type I error
#data_plot=bind_rows(summary_type1error_power[3],summary_type1error_power[4])
data_plot= bind_rows(summary_delta[c(1,3)],summary_m[c(1,3),],summary_f[c(1,3),])
data_plot$method = c(rep(c("PGS-TRI (DE+IDE)","Logistic Regression + PGS(MF)"),2),rep("Logistic Regression + PGS(MF)",4))
data_plot$method=factor(data_plot$method,levels=c("PGS-TRI (DE+IDE)","Logistic Regression + PGS(MF)"))

data_plot$truedelta=c(rep("\u03B2(G)=0,\u03B4-IDE=0",2),rep("\u03B2(G)=0.4,\u03B4-IDE=0",2),
                      rep(c("\u03B2(G)=0,\u03B2(M)=0,\u03B2(F) = 0","\u03B2(G)=0.4,\u03B2(M)=0,\u03B2(F)=0"),2))

data_plot$truedelta=factor(data_plot$truedelta,levels=unique(data_plot$truedelta))

colnames(data_plot)[7]="Method"
data_plot$type = c(rep("\u03B4-IDE",4),rep("\u03B2(Mother)",2),rep("\u03B2(Father)",2))
data_plot$type=factor(data_plot$type,levels = c("\u03B2(Mother)","\u03B2(Father)","\u03B4-IDE"))
data_plot$effect = "\u03B4-IDE"

p2 <- ggplot(data=data_plot, aes(x=truedelta, y=type1error.power , fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width = 0.9,linewidth = 0.1)+
  scale_fill_manual(values = c(pal_nejm("default")(5)[c(3)],"Light Sky Blue")) +
  scale_x_discrete(NULL) +
  scale_y_continuous('Type I Error for \u03B4-IDE',breaks = c(0.05,0.1),limits = c(0,0.1))+
  geom_hline(yintercept=0.05, col = "coral",linetype = "dashed",size=0.5)  + 
  facet_grid(effect~type,scales = "free") + theme_bw()+      
  theme(plot.title = element_text(color = "black", size = 7, face = "bold")) +                                                          
  My_Theme+theme(legend.position="") 


#####################################################################
#####################################################################
#####################################################################
# Panel C


#plot bias and sd

data_plot=bind_rows(summary_beta[c(1,3)]) 
data_plot = data_plot[which(data_plot$method %in% c("logistic","fam","logistic_parent_adj","fam_parent_adj")),]

data_plot$method = rep(c("PGS-TRI (DE only)","Logistic Regression","PGS-TRI (DE+IDE)","Logistic Regression + PGS(MF)"),2)
data_plot$method=factor(data_plot$method,levels=c("PGS-TRI (DE only)","PGS-TRI (DE+IDE)","Logistic Regression","Logistic Regression + PGS(MF)"))

colnames(data_plot)[7]="Method"
data_plot$effect = "DE"


data_plot$truebeta=c(rep("\u03B2(G)=0,\u03B4-IDE=0",4),rep("\u03B2(G)=0.4,\u03B4-IDE=0",4))
data_plot$truebeta=factor(data_plot$truebeta,levels=c("\u03B2(G)=0,\u03B4-IDE=0","\u03B2(G)=0.4,\u03B4-IDE=0"))

data_plot$effect = "DE"

p3 <- ggplot(data_plot, aes(x = truebeta,y=bias,
                            ymin = bias - sd,
                            ymax = bias + sd,colour = Method , fill = Method)) +
  scale_x_discrete(NULL) +
  scale_y_continuous('Bias +- SD \n',lim = c(-0.2, 0.2)
  ) +
  # coord_cartesian(ylim = c(-0.5, 0.5)) +
  # theme_classic() +
  # theme(panel.border = element_rect(fill=NA)) +
  geom_errorbar(position = position_dodge(width=0.6),size=1,width=0.4, aes(color=Method))+ scale_color_npg()  +
  geom_hline(yintercept = 0, colour = "grey", linetype='dashed', size=0.5) +
  scale_color_manual(values = c(pal_nejm("default")(5)[c(1,3)],"Royal Blue","Light Sky Blue",pal_nejm("default")(5)[c(4)])) +
  geom_point(shape = 19, size = 2,position = position_dodge(width=0.6)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  facet_grid(~effect)+ theme_bw()+      
  theme(plot.title = element_text(color = "black", size = 7, face = "bold")) +                                                          
  My_Theme+theme(legend.position="") 

# Panel D


#plot bias and sd

load("./family/simulation/snipar/08192025/results_summary_delta_indirect.rda")
load("./family/simulation/snipar/08192025/results_summary_mother.rda")
load("./family/simulation/snipar/08192025/results_summary_father.rda")

data_plot= bind_rows(summary_delta[c(1,3)],summary_m[c(1,3),],summary_f[c(1,3),])
data_plot$method = c(rep(c("PGS-TRI (DE+IDE)","Logistic Regression + PGS(MF)"),2),rep("Logistic Regression + PGS(MF)",4))
data_plot$method=factor(data_plot$method,levels=c("PGS-TRI (DE+IDE)","Logistic Regression + PGS(MF)"))

data_plot$truedelta=c(rep("\u03B2(G)=0,\u03B4-IDE=0",2),rep("\u03B2(G)=0.4,\u03B4-IDE=0",2),
                      rep(c("\u03B2(G)=0,\u03B2(M)=0,\u03B2(F) = 0","\u03B2(G)=0.4,\u03B2(M)=0,\u03B2(F)=0"),2))


colnames(data_plot)[7]="Method"
data_plot$type = c(rep("\u03B4-IDE",4),rep("\u03B2(Mother)",2),rep("\u03B2(Father)",2))
data_plot$type=factor(data_plot$type,levels = c("\u03B2(Mother)","\u03B2(Father)","\u03B4-IDE"))
data_plot$effect = "\u03B4-IDE"


p4 <- ggplot(data_plot, aes(x = truedelta, bias,
                            ymin = bias - sd,
                            ymax = bias + sd,colour = Method , fill = Method)) +
  scale_x_discrete(NULL) +
  scale_y_continuous('Bias +- SD \n',lim = c(-0.2, 0.2)
  ) +
  geom_errorbar(position = position_dodge(width=0.6),size=1,width=0.2, aes(color=Method))+ 
  geom_hline(yintercept = 0, colour = "grey", linetype='dashed', size=0.5) +
  scale_color_manual(values = c(pal_nejm("default")(5)[c(3)],"Light Sky Blue")) +
  geom_point(shape = 19, size = 2,position = position_dodge(width=0.6)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  facet_grid(effect~type,scales = "free")+ theme_bw()+      
  theme(plot.title = element_text(color = "black", size = 7, face = "bold")) +                                                          
  My_Theme+theme(legend.position="") 




#####################################################################
#####################################################################
#####################################################################
#Legends
#legend_b <- cowplot::get_legend(p1 + theme(legend.position="bottom"))

p_tmp <- p1+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size = 7),legend.title = element_text(face = "bold", size = 7),legend.key.size = unit(3, 'mm'))

legend_b <- cowplot::get_plot_component(p_tmp+ theme(legend.position="top"), 'guide-box-top', return_all = TRUE)


#####################################################################
#####################################################################
#####################################################################


#Finally, put the plots together


p1 <- p1+ theme(legend.position="none")
p2 <- p2+ theme(legend.position="none")
p3 <- p3+ theme(legend.position="none")
p4 <- p4+ theme(legend.position="none")

pg1 <- ggplotGrob(p1)
pg2 <- ggplotGrob(p2)
pg3 <- ggplotGrob(p3)
pg4 <- ggplotGrob(p4)

pg2$heights = pg1$heights
pg3$heights = pg4$heights

p <- ggarrange(ggarrange(p1, p2, 
                         ncol = 2, labels = c("a", "b"),font.label = list(size = 7),
                         widths = c(0.4,0.6),align = "h") , 
               ggarrange(pg3, pg4,
                         ncol = 2, labels = c("c", "d"),font.label = list(size = 7),
                         widths = c(0.4,0.6)) ,
               nrow = 2, 
               labels = c(NA, NA),font.label = list(size = 7),
               heights = c(0.5,0.5))

p <- cowplot::plot_grid( legend_b, p,  ncol = 1, rel_heights = c(.05, 1))

ggsave(filename=paste0("Extended_Figure5.pdf"),
       plot=p, 
       path="./family/manuscript/draft4/plots/extended/",
       width=180, height=160, units="mm", dpi=320,device = cairo_pdf)

ggsave(filename=paste0("Extended_Figure5.png"),
       plot=p, 
       path="./family/manuscript/draft4/plots/extended/",
       width=180, height=160, units="mm", dpi=320)

