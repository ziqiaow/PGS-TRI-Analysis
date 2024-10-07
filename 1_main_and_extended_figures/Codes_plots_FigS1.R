

#----------------------------------------------------------
#Figure S1
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
# Panel A

#95% Coverage Probability for PGS main effect
load("./family/simulation/04262023/06112024/results_summary.rda")
data_plot=bind_rows(summary_ci)

data_plot$truebeta=c(rep("\u03B2(PGS)=0",12),rep("\u03B2(PGS)=0",12),rep("\u03B2(PGS)=0.4",12),
                     rep("\u03B2(PGS)=0.4",12))
data_plot$cor=c(rep("Cor(\u03b1,\u03bc) = 0",12),rep("Cor(\u03b1,\u03bc) = 0.25",12),rep("Cor(\u03b1,\u03bc) = 0",12),rep("Cor(\u03b1,\u03bc) = 0.25",12))

data_plot$cor=factor(data_plot$cor,levels=unique(data_plot$cor))
data_plot$truebeta=factor(data_plot$truebeta,levels=c("\u03B2(PGS)=0","\u03B2(PGS)=0.4"))
colnames(data_plot)[2]="Method"
colnames(data_plot)[3]="Family_Size"
data_plot$Family_Size=as.character(data_plot$Family_Size)
data_plot$Family_Size=factor(data_plot$Family_Size,level=c("200","500","1000","2000"))
data_plot$Method[which(data_plot$Method=="fam")]="PGS-TRI"
data_plot$Method[which(data_plot$Method=="logistic")]="Logistic Regression"
data_plot$Method[which(data_plot$Method=="ptdt")]="pTDT"
data_plot$Method=factor(data_plot$Method,level=c("PGS-TRI","Logistic Regression","pTDT"))
data_plot$beta_prs=data_plot$cp_beta_prs*100

p1 <- ggplot(data=data_plot, aes(x=Family_Size, y=beta_prs , fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width = 0.7)+
  scale_fill_manual(values = c(pal_nejm("default")(5)[c(1)],"Royal Blue",pal_nejm("default")(5)[c(4)])) +
  scale_x_discrete('# of Families') +
  scale_y_continuous('95% Coverage Probability (%)',breaks = c(0,50,95),limits = c(0,100))+
  geom_hline(yintercept=95, col = "coral",linetype = "dashed",size=2)  + 
  facet_grid(cor~truebeta)+ #,labeller = labeller(beta_pgs = supp.labs)
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+ labs(title="Direct PGS Effect") +
  theme(plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme_bw()+    theme(plot.title = element_text(color = "black", size = 12, face = "bold")) +                                                            # Change font size
  theme(strip.text.x = element_text(face = "bold", size = 12),strip.text.y = element_text(face = "bold", size = 12))+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size=12),legend.title = element_text(face = "bold", size=12))+theme(axis.title.x = element_text(face="bold",size=12, color="black"), axis.text.x = element_text(face="bold",size=12, color="black"), axis.title.y = element_text(face="bold",size=12, color="black"), axis.text.y = element_text(face="bold", size=10, color="black"))


#####################################################################
#####################################################################
#####################################################################
# Panel B


#95% Coverage Probability for IDE

load("./family/simulation/nurture/redo/05212024/results_summary.rda")
data_plot=bind_rows(summary_delta)
data_plot$cor=c(rep("Cor(\u03b1,\u03bc) = 0.25",24),rep("Cor(\u03b1,\u03bc) = 0",24))
data_plot$truedelta=rep(c(rep("\u03B4-IDE = 0",8),rep("\u03B4-IDE = 0.2",8),rep("\u03B4-IDE = -0.2",8)),2)

data_plot$cor=factor(data_plot$cor,levels=c("Cor(\u03b1,\u03bc) = 0","Cor(\u03b1,\u03bc) = 0.25"))
data_plot$truedelta=factor(data_plot$truedelta,levels=unique(data_plot$truedelta))
colnames(data_plot)[7]="Method"
colnames(data_plot)[8]="Family_Size"
data_plot$Family_Size=as.character(data_plot$Family_Size)
data_plot$Family_Size=factor(data_plot$Family_Size,level=c("200","500","1000","2000"))
data_plot$Method[which(data_plot$Method=="fam")]="PGS-TRI"
data_plot$Method[which(data_plot$Method=="logistic")]="Logistic Regression"
data_plot$Method=factor(data_plot$Method,level=c("PGS-TRI","Logistic Regression"))

data_plot$coverage_prob=data_plot$coverage_prob*100



p2 <- ggplot(data=data_plot, aes(x=Family_Size, y=coverage_prob , fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width = 0.9)+
  scale_fill_manual(values = c(pal_nejm("default")(5)[c(1)],"Royal Blue")) +
  scale_x_discrete('# of Families') +
  scale_y_continuous('95% Coverage Probability (%)',breaks = c(0,50,95),limits = c(0,100))+
  geom_hline(yintercept=95, col = "coral",linetype = "dashed",size=2)  + 
  facet_grid(cor~truedelta)+ #,labeller = labeller(beta_pgs = supp.labs)
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+   labs(title="Parental Differential Indirect PGS Effect") +
  
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme_bw()+  theme(plot.title = element_text(color = "black", size = 12, face = "bold")) +                                                              # Change font size
  theme(strip.text.x = element_text(face = "bold", size = 12),strip.text.y = element_text(face = "bold", size = 12))+theme(legend.position="top")+  theme(legend.text = element_text(face = "bold", size=12),legend.title = element_text(face = "bold", size=12))+theme(axis.title.x = element_text(face="bold",size=12, color="black"), axis.text.x = element_text(face="bold",size=12, color="black"), axis.title.y = element_text(face="bold",size=12, color="black"), axis.text.y = element_text(face="bold", size=10, color="black"))


#####################################################################
#####################################################################
#####################################################################
# Panel C and D


#95% Coverage Probability for PGSxE1 and PGSxE2

load("./family/simulation_gxe/07232024/cor0/results_summary.rda")

data_plot=bind_rows(summary_ci[c(1,3,5,7)])

data_plot$cor=c(rep("PGS Indep of E",12*2),rep("Cor(PGS,E) \u2260 0",12*2))
data_plot$strat = c(rep("Cor(\u03b1,\u03bc) = 0",12),rep("Cor(\u03b1,\u03bc)=Cor(\u03b1,\u03bc\u03b3)=0.25",12),rep("Cor(\u03b1,\u03bc) = 0",12),rep("Cor(\u03b1,\u03bc)=Cor(\u03b1,\u03bc\u03b3)=0.25",12))

data_plot$cor=factor(data_plot$cor,levels=unique(data_plot$cor))
data_plot$strat=factor(data_plot$strat,levels=unique(data_plot$strat))

colnames(data_plot)[4]="Method"
colnames(data_plot)[5]="Family_Size"
data_plot$Family_Size=as.character(data_plot$Family_Size)
data_plot$Family_Size=factor(data_plot$Family_Size,level=c("200","500","1000","2000"))
data_plot$Method[which(data_plot$Method=="fam")]="PGS-TRI"
data_plot$Method[which(data_plot$Method=="logistic")]="Logistic Regression"
data_plot$Method[which(data_plot$Method=="caseonly")]="Case-Only"
data_plot$Method=factor(data_plot$Method,level=c("PGS-TRI","Logistic Regression","Case-Only"))
data_plot$beta_prs=data_plot$mse_beta_prs*100
data_plot$beta_prs_e1=data_plot$beta_prs_e1*100
data_plot$beta_prs_e2=data_plot$beta_prs_e2*100


p3 <- ggplot(data=data_plot, aes(x=Family_Size, y=beta_prs_e1 , fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width = 0.9)+
  scale_fill_manual(values = c(pal_nejm("default")(5)[c(1)],"Royal Blue",pal_jco("default")(3)[2]))+
  scale_x_discrete('# of Families') +
  scale_y_continuous('95% Coverage Probability (%)',breaks = c(0,50,95),limits = c(0,100))+
  geom_hline(yintercept=95, col = "coral",linetype = "dashed",size=2)  + 
  facet_grid(strat~cor)+ #,labeller = labeller(beta_pgs = supp.labs)
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+  labs(title="PGS x E1 Interaction") +
  theme(plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme_bw()+           theme(plot.title = element_text(color = "black", size = 12, face = "bold")) +                                                     # Change font size
  theme(strip.text.x = element_text(face = "bold", size = 12),strip.text.y = element_text(face = "bold", size = 12))+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size=12),legend.title = element_text(face = "bold", size=12))+theme(axis.title.x = element_text(face="bold",size=12, color="black"), axis.text.x = element_text(face="bold",size=12, color="black"), axis.title.y = element_text(face="bold",size=12, color="black"), axis.text.y = element_text(face="bold", size=10, color="black"))


p4 <- ggplot(data=data_plot, aes(x=Family_Size, y=beta_prs_e2 , fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width = 0.9)+
  scale_fill_manual(values = c(pal_nejm("default")(5)[c(1)],"Royal Blue",pal_jco("default")(3)[2]))+
  scale_x_discrete('# of Families') +
  scale_y_continuous('95% Coverage Probability (%)',breaks = c(0,50,95),limits = c(0,100))+
  geom_hline(yintercept=95, col = "coral",linetype = "dashed",size=2)  + 
  facet_grid(strat~cor)+ #,labeller = labeller(beta_pgs = supp.labs)
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+  labs(title="PGS x E2 Interaction") +
  theme(plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme_bw()+      theme(plot.title = element_text(color = "black", size = 12, face = "bold")) +                                                          # Change font size
  theme(strip.text.x = element_text(face = "bold", size = 12),strip.text.y = element_text(face = "bold", size = 12))+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size=12),legend.title = element_text(face = "bold", size=12))+theme(axis.title.x = element_text(face="bold",size=12, color="black"), axis.text.x = element_text(face="bold",size=12, color="black"), axis.title.y = element_text(face="bold",size=12, color="black"), axis.text.y = element_text(face="bold", size=10, color="black"))




#####################################################################
#####################################################################
#####################################################################
# legends
data_plot$Method = as.character(data_plot$Method)
data_plot$Method[7] = "pTDT"
data_plot$Method=factor(data_plot$Method,level=c("PGS-TRI","Logistic Regression","pTDT","Case-Only"))


p_tmp <-  ggplot(data=data_plot, aes(x=Family_Size, y=beta_prs_e2 , fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width = 0.9)+
  scale_fill_manual(values = c(pal_nejm("default")(5)[1],"Royal Blue",pal_nejm("default")(5)[c(4)],pal_jco("default")(3)[2]))+
  scale_x_discrete('# of Families') +
  scale_y_continuous(NULL,breaks = c(0,50,95),limits = c(0,100))+
  geom_hline(yintercept=95, col = "coral",linetype = "dashed",size=2)  + 
  facet_grid(strat~cor)+ #,labeller = labeller(beta_pgs = supp.labs)
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+  labs(title="PGS x E2 Interaction") +
  theme(plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme_bw()+                                                                # Change font size
  theme(strip.text.x = element_text(face = "bold", size = 12),strip.text.y = element_text(face = "bold", size = 12))+theme(legend.position="right")+  theme(legend.text = element_text(face = "bold", size=12),legend.title = element_text(face = "bold", size=12))+theme(axis.title.x = element_text(face="bold",size=12, color="black"), axis.text.x = element_text(face="bold",size=12, color="black"), axis.title.y = element_text(face="bold",size=12, color="black"), axis.text.y = element_text(face="bold", size=10, color="black"))

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



p <- ggarrange(ggarrange(p1,p2,
                         ncol = 2, labels = c("a", "b"),
                         widths = c(0.5,0.5)) , 
               ggarrange(p3,p4,
                         ncol = 2, labels = c("c","d"),
                         widths = c(0.5,0.5)) ,
               nrow = 2, 
               labels = c(NA, NA),
               heights = c(0.5,0.5))


p <- plot_grid( legend_b,p,  ncol = 1, rel_heights = c( .05,1))

ggsave(filename=paste0("Extended_Figure1_09262024.pdf"),
       plot=p, 
       path="./family/manuscript/plots/final/",
       width=360, height=350, units="mm", dpi=320,device = cairo_pdf)

ggsave(filename=paste0("Extended_Figure1_09262024.png"),
       plot=p, 
       path="./family/manuscript/plots/final/",
       width=360, height=350, units="mm", dpi=320)
