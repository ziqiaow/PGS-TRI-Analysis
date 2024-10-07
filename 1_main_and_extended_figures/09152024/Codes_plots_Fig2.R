

#----------------------------------------------------------
#Figure 2
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



#####################################################################
#####################################################################
#####################################################################
# Panel A

#type I error

#Load data for DE
load("./family/simulation/04262023/06112024/results_summary.rda")

#Load data for IDE
load("./family/simulation/nurture/redo/05212024/results_summary.rda")

data_plot=bind_rows(summary_type1error_power[c(1,2)])
data_plot_IDE = bind_rows(summary_delta[c(4,1)])
data_plot_IDE = data_plot_IDE[,c(4,7,8)]
colnames(data_plot_IDE) = colnames(data_plot)
data_plot_IDE$measure = "\u03B4-IDE"
data_plot$measure = "DE"
data_plot = rbind(data_plot,data_plot_IDE)


data_plot = data_plot[which(data_plot$samplesize == "200"),]
data_plot$cor=c(rep("Cor(\u03b1,\u03bc) = 0",3),rep("Cor(\u03b1,\u03bc) = 0.25",3),rep("Cor(\u03b1,\u03bc) = 0",2),rep("Cor(\u03b1,\u03bc) = 0.25",2))



data_plot$cor=factor(data_plot$cor,levels=unique(data_plot$cor))

colnames(data_plot)[2]="Method"
colnames(data_plot)[3]="Family_Size"
data_plot$Method[which(data_plot$Method=="fam")]="PGS-TRI"
data_plot$Method[which(data_plot$Method=="logistic")]="Logistic Regression"
data_plot$Method[which(data_plot$Method=="ptdt")]="pTDT"
data_plot$Method=factor(data_plot$Method,level=c("PGS-TRI","Logistic Regression","pTDT"))

p1 <- ggplot(data=data_plot, aes(x=Method, y=type1error_beta_prs , fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width = 0.9)+
  geom_hline(aes(yintercept=0.15), alpha = 0)+
  scale_fill_manual(values = c(pal_nejm("default")(5)[c(1)],"Royal Blue",pal_nejm("default")(5)[c(4)])) +
  scale_x_discrete('') +ggtitle("Direct and Indirect PGS Effect")+
  scale_y_continuous('Type I Error',breaks = c(0.05,0.5,1),limits = c(0,1))+
  geom_hline(yintercept=0.05, col = "coral",linetype = "dashed",size=2)  + 
  facet_grid(cor~measure,scales="free")+ 
  theme(axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 14))+ 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme_bw()+      theme(plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5)) +                                                           # Change font size
  theme(strip.text.x = element_text(face = "bold", size = 14),strip.text.y = element_text(face = "bold", size = 14))+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size = 14),legend.title = element_text(face = "bold", size = 14))+theme(axis.title.x = element_text(face="bold",size = 14, color="black"), axis.text.x = element_text(face="bold",size = 14, color="black",angle = 45, vjust = 0.5, hjust = 0.5), axis.title.y = element_text(face="bold",size = 14, color="black"), axis.text.y = element_text(face="bold", size = 14, color="black"))







#####################################################################
#####################################################################
#####################################################################
# Panel C


#plot bias and sd

#Load data for DE
load("./family/simulation/04262023/06112024/results_summary.rda")

#Load data for IDE
load("./family/simulation/nurture/redo/05212024/results_summary.rda")

data_plot_IDE = bind_rows(summary_delta[c(1,2,4,5)])
data_plot_IDE = data_plot_IDE[,c(1,2,7,8,10)]
data_plot_IDE$cor=c(rep("Cor(\u03b1,\u03bc) = 0.25",16),rep("Cor(\u03b1,\u03bc) = 0",16))
data_plot_IDE = data_plot_IDE[which(data_plot_IDE$family_size == "1000"),]
data_plot_IDE$measure = "\u03B4-IDE"
data_plot_IDE$true_delta=c(rep("\u03B4-IDE=0",2),rep("\u03B4-IDE=0.2",2),rep("\u03B4-IDE=0",2),rep("\u03B4-IDE=0.2",2))


data_plot=bind_rows(summary_bias)
data_plot$truebeta=c(rep("\u03B2(G)=0",24),rep("\u03B2(G)=0.4",24))
data_plot$cor=c(rep("Cor(\u03b1,\u03bc) = 0",12),rep("Cor(\u03b1,\u03bc) = 0.25",12),rep("Cor(\u03b1,\u03bc) = 0",12),rep("Cor(\u03b1,\u03bc) = 0.25",12))
colnames(data_plot)[2]="Method"
colnames(data_plot)[3]="Family_Size"
data_plot = data_plot[-which(data_plot$Method=="ptdt"),]
data_plot = data_plot[which(data_plot$Family_Size == "1000"),]

#add SD
data_plot2=bind_rows(summary_empirical_sd)
data_plot2= data_plot2[-which(data_plot2$method == "ptdt"),]
data_plot2 = data_plot2[which(data_plot2$samplesize == "1000"),]
data_plot=cbind(data_plot,data_plot2[,1])
colnames(data_plot)[6]=c("sd1")
data_plot = data_plot[,c(1,6,2,3,4,5)]
data_plot$measure = "DE"
colnames(data_plot_IDE) = colnames(data_plot)
data_plot = rbind(data_plot,data_plot_IDE)

data_plot$cor=factor(data_plot$cor,levels=unique(data_plot$cor))
data_plot$truebeta = as.character(data_plot$truebeta)
data_plot$truebeta=factor(data_plot$truebeta,levels=c("\u03B2(G)=0","\u03B2(G)=0.4","\u03B4-IDE=0","\u03B4-IDE=0.2"))

data_plot$Method[which(data_plot$Method=="fam")]="PGS-TRI"
data_plot$Method[which(data_plot$Method=="logistic")]="Logistic Regression"
data_plot$Method=factor(data_plot$Method,level=c("PGS-TRI","Logistic Regression"))



p3 <- ggplot(data_plot, aes(x = truebeta, bias_beta_prs,
                      ymin = bias_beta_prs - sd1,
                      ymax = bias_beta_prs + sd1,colour = Method , fill = Method)) +
  scale_x_discrete(NULL) +ggtitle("Direct and Indirect PGS Effect")+
  scale_y_continuous('Bias +- SD',lim = c(-0.3, 0.3)
  ) +
  geom_errorbar(position = position_dodge(width=0.6),size=1,width=0.4, aes(color=Method))+ scale_color_npg()  +
  geom_hline(yintercept=0, col = "grey")  +
  scale_color_manual(values = c(pal_nejm("default")(5)[c(1)],"Royal Blue",pal_nejm("default")(5)[c(4)])) +
  geom_point(shape = 19, size = 4,position = position_dodge(width=0.6)) +
  # theme(legend.position="none") +
  theme(axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  facet_grid(cor~measure,scales="free_x")+ theme_bw()+ theme(plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5)) +                                                         # Change font size
  theme(strip.text.x = element_text(face = "bold", size = 14),strip.text.y = element_text(face = "bold", size = 14))+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size = 14),legend.title = element_text(face = "bold", size = 14))+theme(axis.title.x = element_text(face="bold",size = 16, color="black"), axis.text.x = element_text(face="bold",size = 14, color="black",angle = 45, vjust = 0.5, hjust = 0.5), axis.title.y = element_text(face="bold",size = 14, color="black"), axis.text.y = element_text(face="bold", size = 14, color="black"))







#####################################################################
#####################################################################
#####################################################################
# Panel B

#type I error for PGSxE interactions
load("./family/simulation_gxe/07232024/cor0/results_summary.rda")
data_plot=bind_rows(summary_type1error_power[c(1,3,5,7)])
#select when trios = 1000
data_plot = data_plot[which(data_plot$samplesize == "1000"),]
data_plot$cor=c(rep("PGS Indep of E",6),rep("Cor(PGS,E) \u2260 0",6))
data_plot$strat = c(rep("Cor(\u03b1,\u03bc) = 0",3),rep("Cor(\u03b1,\u03bc)=Cor(\u03b1,\u03bc\u03b3)=0.25",3),rep("Cor(\u03b1,\u03bc) = 0",3),rep("Cor(\u03b1,\u03bc)=Cor(\u03b1,\u03bc\u03b3)=0.25",3))

data_plot$cor=factor(data_plot$cor,levels=unique(data_plot$cor))
data_plot$strat=factor(data_plot$strat,levels=unique(data_plot$strat))


colnames(data_plot)[4]="Method"
data_plot$Method[which(data_plot$Method=="fam")]="PGS-TRI"
data_plot$Method[which(data_plot$Method=="logistic")]="Logistic Regression"
data_plot$Method[which(data_plot$Method=="caseonly")]="Case-Only"
data_plot$Method=factor(data_plot$Method,level=c("PGS-TRI","Logistic Regression","Case-Only"))


p2 <- ggplot(data=data_plot, aes(x=Method, y=beta_prs_e1 , fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width = 0.9)+
  scale_fill_manual(values = c(pal_nejm("default")(5)[1],"Royal Blue",pal_jco("default")(3)[2]))+
  scale_x_discrete(NULL) +ggtitle("PGS x E Interactions")+
  scale_y_continuous('Type I Error',breaks = c(0.05,0.5,1))+
  geom_hline(yintercept=0.05, col = "coral",linetype = "dashed",size=2)  + 
  facet_grid(strat~cor)+ #,labeller = labeller(beta_pgs = supp.labs)
  theme(axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 14))+  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme_bw()+         theme(plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5)) +                                                         # Change font size
theme(strip.text.x = element_text(face = "bold", size = 14),strip.text.y = element_text(face = "bold", size = 11))+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size = 14),legend.title = element_text(face = "bold", size = 14))+theme(axis.title.x = element_text(face="bold",size = 14, color="black"), axis.text.x = element_text(face="bold",size = 14, color="black",angle = 45, vjust = 0.5, hjust = 0.5), axis.title.y = element_text(face="bold",size = 14, color="black"), axis.text.y = element_text(face="bold", size = 14, color="black"))




#####################################################################
#####################################################################
#####################################################################
# Panel D


#plot bias and sd for PGSxE
load("./family/simulation_gxe/07232024/cor0/results_summary.rda")

data_plot=bind_rows(summary_bias)
data_plot$truebeta=c(rep("\u03B2(GxE)=0",12),rep("\u03B2(GxE)=0.3",12),rep("\u03B2(GxE)=0",12),rep("\u03B2(GxE)=0.3",12),
                     rep("\u03B2(GxE)=0",12),rep("\u03B2(GxE)=0.3",12),rep("\u03B2(GxE)=0",12),rep("\u03B2(GxE)=0.3",12))
data_plot = data_plot[which(data_plot$samplesize == "1000"),]
data_plot$cor=c(rep("PGS Indep of E",12),rep("Cor(PGS,E) \u2260 0",12))
data_plot$strat = c(rep("Cor(\u03b1,\u03bc) = 0",6),rep("Cor(\u03b1,\u03bc)=Cor(\u03b1,\u03bc\u03b3)=0.25",6),rep("Cor(\u03b1,\u03bc) = 0",6),rep("Cor(\u03b1,\u03bc)=Cor(\u03b1,\u03bc\u03b3)=0.25",6))

data_plot$cor=factor(data_plot$cor,levels=unique(data_plot$cor))
data_plot$strat=factor(data_plot$strat,levels=unique(data_plot$strat))
colnames(data_plot)[4]="Method"
data_plot$Method[which(data_plot$Method=="fam")]="PGS-TRI"
data_plot$Method[which(data_plot$Method=="logistic")]="Logistic Regression"
data_plot$Method[which(data_plot$Method=="caseonly")]="Case-Only"
data_plot$Method=factor(data_plot$Method,level=c("PGS-TRI","Logistic Regression","Case-Only"))

data_plot2=bind_rows(summary_empirical_sd)
data_plot2 = data_plot2[which(data_plot2$samplesize == "1000"),]
data_plot2=data_plot2[,-c(4,5)]
colnames(data_plot2)=c("sd1","sd2","sd3")

data_plot=cbind(data_plot,data_plot2)

p4 <- ggplot(data_plot, aes(x = truebeta, beta_prs_e1,
                      ymin = beta_prs_e1 - sd2,
                      ymax = beta_prs_e1 + sd2,colour = Method , fill = Method)) +
  scale_x_discrete(NULL) +ggtitle("PGS x E Interactions")+
  scale_y_continuous('Bias +- SD',lim = c(-0.3, 0.3)
  ) +
  geom_errorbar(position = position_dodge(width=0.6),size=1,width=0.4, aes(color=Method))+ scale_color_npg()  + 
  geom_hline(yintercept=0, col = "grey")  + 
  geom_point(shape = 19, size = 4,position = position_dodge(width=0.6)) +
  theme(axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 14)) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  facet_grid(strat~cor)+ theme_bw()+   theme(plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5)) +scale_color_manual(values = c(pal_nejm("default")(5)[1],"Royal Blue",pal_jco("default")(3)[2]))+                                                          # Change font size
  theme(strip.text.x = element_text(face = "bold", size = 14),strip.text.y = element_text(face = "bold", size = 11))+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size = 14),legend.title = element_text(face = "bold", size = 14))+theme(axis.title.x = element_text(face="bold",size = 16, color="black"), axis.text.x = element_text(face="bold",size = 14, color="black",angle = 45, vjust = 0.5, hjust = 0.5), axis.title.y = element_text(face="bold",size = 14, color="black"), axis.text.y = element_text(face="bold", size = 14, color="black"))





#####################################################################
#####################################################################
#####################################################################
# legends
data_plot$Method = as.character(data_plot$Method)
data_plot$Method[which(data_plot$Method=="fam")]="PGS-TRI"
data_plot$Method[which(data_plot$Method=="logistic")]="Logistic Regression"
data_plot$Method[which(data_plot$Method=="caseonly")]="Case-Only"
data_plot$Method[6] = "pTDT"
data_plot$Method=factor(data_plot$Method,level=c("PGS-TRI","Logistic Regression","pTDT","Case-Only"))


p_tmp <- ggplot(data_plot, aes(x = truebeta, beta_prs_e1,
                            ymin = beta_prs_e1 - sd2,
                            ymax = beta_prs_e1 + sd2,colour = Method , fill = Method)) +
  scale_x_discrete('# of Families') +
  scale_y_continuous('Bias +- SD \n',lim = c(-0.3, 0.3))+
  geom_errorbar(position = position_dodge(width=0.6),size=2,width=0.4, aes(color=Method))+ 
  geom_hline(yintercept=0, col = "grey")  + 
  geom_point(shape = 19, size = 5,position = position_dodge(width=0.6)) +
  theme(axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 14)) + 
  theme(plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  facet_grid(strat~cor)+ theme_bw()+ scale_color_manual(values = c(pal_nejm("default")(5)[1],"Royal Blue",pal_nejm("default")(5)[c(4)],pal_jco("default")(3)[2]))+    
# Change font size
  theme(strip.text.x = element_text(face = "bold", size = 14),strip.text.y = element_text(face = "bold", size = 14))+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size = 14),legend.title = element_text(face = "bold", size = 14))+theme(axis.title.x = element_text(face="bold",size = 14, color="black"), axis.text.x = element_text(face="bold",size = 14, color="black"), axis.title.y = element_text(face="bold",size = 14, color="black"), axis.text.y = element_text(face="bold", size = 14, color="black"))


p_tmp <- p_tmp+  theme(strip.text.x = element_text(face = "bold", size = 18),strip.text.y = element_text(face = "bold", size = 18))+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size = 18),legend.title = element_text(face = "bold", size = 18))+theme(axis.title.x = element_text(face="bold",size = 18, color="black"), axis.text.x = element_text(face="bold",size = 18, color="black",angle = 45, vjust = 0.5, hjust = 0.5), axis.title.y = element_text(face="bold",size = 18, color="black"), axis.text.y = element_text(face="bold", size = 18, color="black"))

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

pg2$heights <- pg1$heights

p <- ggarrange(ggarrange(pg1, pg2, 
                         ncol = 2, labels = c("a", "b"),
                         widths = c(0.5,0.5)) , 
               ggarrange(p3, p4,
                         ncol = 2, labels = c("c", "d"),
                         widths = c(0.5,0.5)) ,
               nrow = 2, 
               labels = c(NA, NA),
               heights = c(0.53,0.47))

p <- cowplot::plot_grid( legend_b, p,  ncol = 1, rel_heights = c(.05, 1))

ggsave(filename=paste0("Figure2_addGxE_4panels_09262024.pdf"),
       plot=p, 
       path="./family/manuscript/plots/final/",
       width=360, height=350, units="mm", dpi=320,device = cairo_pdf)

ggsave(filename=paste0("Figure2_addGxE_4panels_09262024.png"),
       plot=p, 
       path="./family/manuscript/plots/final/",
       width=360, height=350, units="mm", dpi=320)
