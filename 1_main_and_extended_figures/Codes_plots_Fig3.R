

#----------------------------------------------------------
#Figure 3
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

#type I error for UKB simulation results for DE
#load data with PS matched by geographical regions only

load("./family/simulation/UKB/final_06122024/PS/PGS_main/results_summary.rda")

#type I error
data_plot=bind_rows(summary_type1error_power[4],summary_type1error_power[5])
data_plot = data_plot[-which(data_plot$method == "pcfather" | data_plot$method == "pcmother"| data_plot$method == "logistic_pcfamo"),]
data_plot$Family_Size = c(rep("#Trios = 1000",5),rep("#Trios = 2000",5))

data_plot$method = rep(c("PGS-TRI","Logistic Regression","pTDT","PC10","PC10 Geo"),2)
data_plot$method=factor(data_plot$method,levels=c("PGS-TRI","Logistic Regression","PC10","PC10 Geo","pTDT"))


colnames(data_plot)[2]="Method"
data_plot$measure = "PS"
data_plot1 = data_plot

#load data with population stratification bias and assortative mating matched by geographical regions and educational attainments
load("./family/simulation/UKB/final_06122024/PS_AM/PGS_main/results_summary.rda")
#type I error
data_plot=bind_rows(summary_type1error_power[4],summary_type1error_power[5])
data_plot = data_plot[-which(data_plot$method == "pcfather" | data_plot$method == "pcmother"| data_plot$method == "logistic_pcfamo"),]
data_plot$Family_Size = c(rep("#Trios = 1000",5),rep("#Trios = 2000",5))

data_plot$method = rep(c("PGS-TRI","Logistic Regression","pTDT","PC10","PC10 Geo"),2)
data_plot$method=factor(data_plot$method,levels=c("PGS-TRI","Logistic Regression","PC10","PC10 Geo","pTDT"))


colnames(data_plot)[2]="Method"
data_plot$measure = "PS + AM"


data_plot = rbind(data_plot1,data_plot)
data_plot = data_plot[which(data_plot$Family_Size == "#Trios = 2000"),]
data_plot$measure = factor(data_plot$measure,levels = unique(data_plot$measure))
data_plot$effect = "DE"
#data_plot$type = "\u03B2(PGS)"

p1 <- ggplot(data=data_plot, aes(x=Method, y=type1error_beta_prs , fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width = 0.9)+
  scale_fill_manual(values = c(pal_nejm("default")(5)[c(1)],"Royal Blue","Cornflower Blue","Light Sky Blue",pal_nejm("default")(5)[c(4)])) +
  scale_x_discrete(NULL) +
  scale_y_continuous('Type I Error',breaks = c(0,0.05,0.1),limits = c(0,0.15))+
  geom_hline(yintercept=0.05, col = "coral",linetype = "dashed",size=2)  + 
  facet_grid(measure~effect)+ 
  theme(axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 14))+  #labs(title="PGS Main Effect") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme_bw()+                                                                # Change font size
  theme(strip.text.x = element_text(face = "bold", size = 14),strip.text.y = element_text(face = "bold", size = 14))+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size = 14),legend.title = element_text(face = "bold", size = 14))+theme(axis.title.x = element_text(face="bold",size = 14, color="black"), axis.text.x = element_text(face="bold",size = 14, color="black",angle = 45, vjust = 0.5, hjust = 0.5), axis.title.y = element_text(face="bold",size = 14, color="black"), axis.text.y = element_text(face="bold", size = 14, color="black"))



#####################################################################
#####################################################################
#####################################################################
# Panel C


#plot bias and sd for UKB simulation results
#load data with population stratification bias matched by geographical regions only

load("./family/simulation/UKB/final_06122024/PS/PGS_main/results_summary.rda")


#plot bias and sd

data_plot=bind_rows(summary_bias[c(1,2,4,5)]) 
data_plot = data_plot[-which(data_plot$method == "pcfather" | data_plot$method == "pcmother"| data_plot$method == "logistic_pcfamo"),]
data_plot$truebeta=c(rep("\u03B2(G) = 0.4",10),rep("\u03B2(G) = 0",10))
data_plot$Family_Size = c(rep("#Trios = 1000",5),rep("#Trios = 2000",5),rep("#Trios = 1000",5),rep("#Trios = 2000",5))

data_plot$truebeta=factor(data_plot$truebeta,levels=c("\u03B2(G) = 0","\u03B2(G) = 0.4"))
data_plot$method = rep(c("PGS-TRI","Logistic Regression","pTDT","PC10","PC10 Geo"),4)
data_plot$method=factor(data_plot$method,levels=c("PGS-TRI","Logistic Regression","PC10","PC10 Geo","pTDT"))



colnames(data_plot)[2]="Method"

data_plot3=bind_rows(summary_empirical_sd[c(1,2,4,5)])
data_plot3 = data_plot3[-which(data_plot3$method == "pcfather" | data_plot3$method == "pcmother"| data_plot3$method == "logistic_pcfamo"),]

data_plot=cbind(data_plot,data_plot3[,1])
colnames(data_plot)[5]=c("sd1")
data_plot$measure = "PS"

data_plot1 = data_plot


#load data with population stratification bias and assortative mating matched by geographical regions and educational attainments
load("./family/simulation/UKB/final_06122024/PS_AM/PGS_main/results_summary.rda")

#plot bias and sd

data_plot=bind_rows(summary_bias[c(1,2,4,5)])
data_plot = data_plot[-which(data_plot$method == "pcfather" | data_plot$method == "pcmother"| data_plot$method == "logistic_pcfamo"),]
data_plot$truebeta=c(rep("\u03B2(G) = 0.4",10),rep("\u03B2(G) = 0",10))
data_plot$Family_Size = c(rep("#Trios = 1000",5),rep("#Trios = 2000",5),rep("#Trios = 1000",5),rep("#Trios = 2000",5))

data_plot$truebeta=factor(data_plot$truebeta,levels=c("\u03B2(G) = 0","\u03B2(G) = 0.4"))
data_plot$method = rep(c("PGS-TRI","Logistic Regression","pTDT","PC10","PC10 Geo"),4)
data_plot$method=factor(data_plot$method,levels=c("PGS-TRI","Logistic Regression","PC10","PC10 Geo","pTDT"))



colnames(data_plot)[2]="Method"

data_plot3=bind_rows(summary_empirical_sd[c(1,2,4,5)]) 
data_plot3 = data_plot3[-which(data_plot3$method == "pcfather" | data_plot3$method == "pcmother"| data_plot3$method == "logistic_pcfamo"),]

data_plot=cbind(data_plot,data_plot3[,1])
colnames(data_plot)[5]=c("sd1")
data_plot$measure = "PS + AM"


data_plot = rbind(data_plot1,data_plot)
data_plot = data_plot[-which(data_plot$Method == "pTDT"),]
data_plot$Family_Size = factor(data_plot$Family_Size)


#only plot family size = 2000
data_plot = data_plot[which(data_plot$Family_Size=="#Trios = 2000"),]
data_plot$measure = factor(data_plot$measure,levels = unique(data_plot$measure))
data_plot$effect = "DE"
p3 <- ggplot(data_plot, aes(x = truebeta, bias_beta_prs,
                            ymin = bias_beta_prs - sd1,
                            ymax = bias_beta_prs + sd1,colour = Method , fill = Method)) +
  scale_x_discrete(NULL) +
  scale_y_continuous('Bias +- SD \n',lim = c(-0.6, 0.6)
  ) +
  geom_errorbar(position = position_dodge(width=0.6),size=1.5,width=0.4, aes(color=Method))+ scale_color_npg()  +
  geom_hline(yintercept=0, col = "grey")  +
  scale_color_manual(values = c(pal_nejm("default")(5)[c(1)],"Royal Blue","Cornflower Blue","Light Sky Blue",pal_nejm("default")(5)[c(4)])) +
  geom_point(shape = 19, size = 4,position = position_dodge(width=0.6)) +
  theme(axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 14)) + #labs(title="PGS Main Effect") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  facet_grid(measure~effect)+ theme_bw()+                                                                # Change font size
  theme(strip.text.x = element_text(face = "bold", size = 14),strip.text.y = element_text(face = "bold", size = 14))+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size = 14),legend.title = element_text(face = "bold", size = 14))+theme(axis.title.x = element_text(face="bold",size = 14, color="black"),  axis.text.x = element_text(face="bold",size = 14, color="black",angle = 45, vjust = 0.5, hjust = 0.5), axis.title.y = element_text(face="bold",size = 14, color="black"), axis.text.y = element_text(face="bold", size = 14, color="black"))



#####################################################################
#####################################################################
#####################################################################
# Panel B

#type I error for UKB simulation results for IDE
#load data with PS matched by geographical regions only

load("./family/simulation/UKB/final_06122024/PS/nurture/results_summary.rda")
load("./family/simulation/UKB/final_06122024/PS/nurture/results_summary_beta_main.rda")

#type I error
data_plot=bind_rows(summary_delta,summary_beta)
data_plot$truedelta=c(rep(c(rep("\u03B4-IDE = 0",4),rep("\u03B4-IDE = 0.2",4),rep("\u03B4-IDE = -0.2",4)),2),
                      rep(c(rep("\u03B4-IDE = 0",4),rep("\u03B4-IDE = 0.2",4),rep("\u03B4-IDE = -0.2",4)),2))
data_plot$type = c(rep("\u03B4-IDE",24),rep("DE = 0",24))
data_plot$truedelta=factor(data_plot$truedelta,levels=unique(data_plot$truedelta))
colnames(data_plot)[7]="Method"
colnames(data_plot)[9]="Family_Size"
data_plot$Family_Size=as.character(data_plot$Family_Size)
data_plot$Family_Size=factor(data_plot$Family_Size,level=c("1000","2000"))
data_plot$Method = c(rep(c("PGS-TRI","Logistic Regression","PC10","PC10 Geo"),12))
data_plot$Method=factor(data_plot$Method,levels=c("PGS-TRI","Logistic Regression","PC10","PC10 Geo"))
data_plot$coverage_prob=data_plot$coverage_prob*100

#Type I error
data_plot1=data_plot[which(data_plot$true_delta==0 | data_plot$true_beta == 0 ),]

data_plot1$measure = "PS"


#load data with population stratification bias and assortative mating matched by geographical regions and educational attainments
load("./family/simulation/UKB/final_06122024/PS_AM/nurture/results_summary.rda")
load("./family/simulation/UKB/final_06122024/PS_AM/nurture/results_summary_beta_main.rda")

#type I error
data_plot=bind_rows(summary_delta,summary_beta)
data_plot$truedelta=c(rep(c(rep("\u03B4-IDE = 0",4),rep("\u03B4-IDE = 0.2",4),rep("\u03B4-IDE = -0.2",4)),2),
                      rep(c(rep("\u03B4-IDE = 0",4),rep("\u03B4-IDE = 0.2",4),rep("\u03B4-IDE = -0.2",4)),2))
data_plot$type = c(rep("\u03B4-IDE",24),rep("DE = 0",24))
data_plot$truedelta=factor(data_plot$truedelta,levels=unique(data_plot$truedelta))
colnames(data_plot)[7]="Method"
colnames(data_plot)[9]="Family_Size"
data_plot$Family_Size=as.character(data_plot$Family_Size)
data_plot$Family_Size=factor(data_plot$Family_Size,level=c("1000","2000"))
data_plot$Method = c(rep(c("PGS-TRI","Logistic Regression","PC10","PC10 Geo"),12))
data_plot$Method=factor(data_plot$Method,levels=c("PGS-TRI","Logistic Regression","PC10","PC10 Geo"))
data_plot$coverage_prob=data_plot$coverage_prob*100

#Type I error
data_plot2=data_plot[which(data_plot$true_delta==0 | data_plot$true_beta == 0 ),]
data_plot2$measure = "PS + AM"


data_plot = rbind(data_plot1,data_plot2)
data_plot$measure = factor(data_plot$measure,levels = unique(data_plot$measure))

data_plot = data_plot[which(data_plot$Family_Size == "2000"),]
data_plot = data_plot[-which(data_plot$type == "DE = 0"),]


p2 <- ggplot(data=data_plot, aes(x=Method, y=type1error.power , fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width = 0.9)+
  scale_fill_manual(values = c(pal_nejm("default")(5)[c(1)],"Royal Blue","Cornflower Blue","Light Sky Blue")) +
  scale_x_discrete(NULL) +
  scale_y_continuous('Type I Error',breaks = c(0,0.05,0.1),limits = c(0,0.15))+
  geom_hline(yintercept=0.05, col = "coral",linetype = "dashed",size=2)  + 
  facet_grid(measure~type,scales = "free")+ 
  theme(axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 14))+  #labs(title="PGS Main Effect") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme_bw()+                                                                # Change font size
  theme(strip.text.x = element_text(face = "bold", size = 14),strip.text.y = element_text(face = "bold", size = 14))+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size = 14),legend.title = element_text(face = "bold", size = 14))+theme(axis.title.x = element_text(face="bold",size = 14, color="black"), axis.text.x = element_text(face="bold",size = 14, color="black",angle = 45, vjust = 0.5, hjust = 0.5), axis.title.y = element_text(face="bold",size = 14, color="black"), axis.text.y = element_text(face="bold", size = 14, color="black"))


#####################################################################
#####################################################################
#####################################################################
# Panel D
#load data with PS matched by geographical regions only

load("./family/simulation/UKB/final_06122024/PS/nurture/results_summary.rda")
load("./family/simulation/UKB/final_06122024/PS/nurture/results_summary_beta_main.rda")
#type I error
data_plot=bind_rows(summary_delta,summary_beta)
data_plot$truedelta=c(rep(c(rep("\u03B4-IDE = 0",4),rep("\u03B4-IDE = 0.2",4),rep("\u03B4-IDE = -0.2",4)),2),
                      rep(c(rep("\u03B4-IDE = 0",4),rep("\u03B4-IDE = 0.2",4),rep("\u03B4-IDE = -0.2",4)),2))
data_plot$type = c(rep("\u03B4-IDE",24),rep("\u03B2(G) = 0",24))
data_plot$truedelta=factor(data_plot$truedelta,levels=unique(data_plot$truedelta))
colnames(data_plot)[7]="Method"
colnames(data_plot)[9]="Family_Size"
data_plot$Family_Size=as.character(data_plot$Family_Size)
data_plot$Family_Size=factor(data_plot$Family_Size,level=c("1000","2000"))
data_plot$Method = c(rep(c("PGS-TRI","Logistic Regression","PC10","PC10 Geo"),12))
data_plot$Method=factor(data_plot$Method,levels=c("PGS-TRI","Logistic Regression","PC10","PC10 Geo"))
data_plot$coverage_prob=data_plot$coverage_prob*100
data_plot1=data_plot
data_plot1$measure = "PS"




#load data with population stratification bias and assortative mating matched by geographical regions and educational attainments
load("./family/simulation/UKB/final_06122024/PS_AM/nurture/results_summary.rda")
load("./family/simulation/UKB/final_06122024/PS_AM/nurture/results_summary_beta_main.rda")
#type I error
data_plot=bind_rows(summary_delta,summary_beta)
data_plot$truedelta=c(rep(c(rep("\u03B4-IDE = 0",4),rep("\u03B4-IDE = 0.2",4),rep("\u03B4-IDE = -0.2",4)),2),
                      rep(c(rep("\u03B4-IDE = 0",4),rep("\u03B4-IDE = 0.2",4),rep("\u03B4-IDE = -0.2",4)),2))
data_plot$type = c(rep("\u03B4-IDE",24),rep("\u03B2(G) = 0",24))
#use unicode downloaded from https://www.unicode.org/ucd/, saved the code charts to C:\Users\ziqia\Desktop\work\cluster\CodeCharts.pdf, this is on page 37
#https://stackoverflow.com/questions/5293715/how-to-use-greek-symbols-in-ggplot2

data_plot$truedelta=factor(data_plot$truedelta,levels=unique(data_plot$truedelta))
colnames(data_plot)[7]="Method"
colnames(data_plot)[9]="Family_Size"
data_plot$Family_Size=as.character(data_plot$Family_Size)
data_plot$Family_Size=factor(data_plot$Family_Size,level=c("1000","2000"))
data_plot$Method = c(rep(c("PGS-TRI","Logistic Regression","PC10","PC10 Geo"),12))
data_plot$Method=factor(data_plot$Method,levels=c("PGS-TRI","Logistic Regression","PC10","PC10 Geo"))
data_plot$coverage_prob=data_plot$coverage_prob*100
data_plot2=data_plot
data_plot2$measure = "PS + AM"



data_plot = rbind(data_plot1,data_plot2)
data_plot$measure = factor(data_plot$measure,levels = unique(data_plot$measure))

data_plot = data_plot[which(data_plot$Family_Size == "2000"),]
data_plot = data_plot[-which(data_plot$type == "\u03B2(G) = 0"),]
#plot bias and sd
p4 <- ggplot(data_plot, aes(x = truedelta, bias,
                      ymin = bias - sd,
                      ymax = bias + sd,colour = Method , fill = Method)) +
  scale_x_discrete('') +
  scale_y_continuous('Bias +- SD \n',lim = c(-0.6, 0.6)
  ) +
  geom_errorbar(position = position_dodge(width=0.6),size=1.5,width=0.4, aes(color=Method))+ 
  geom_hline(yintercept=0, col = "grey")  + 
  geom_point(shape = 19, size = 4,position = position_dodge(width=0.6)) +
  scale_color_manual(values =  c(pal_nejm("default")(5)[c(1)],"Royal Blue","Cornflower Blue","Light Sky Blue",pal_nejm("default")(5)[c(4)]))+
  theme(axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 14))+ #+ labs(title="Simulation Investigation of the PGS Main Effect in Case-Parent Trio Study") +
  theme(plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  facet_grid(measure~type,scales = "free")+ theme_bw()+                                                                # Change font size
  theme(strip.text.x = element_text(face = "bold", size = 14),strip.text.y = element_text(face = "bold", size = 14))+theme(legend.position="top")+  theme(legend.text = element_text(face = "bold", size = 14),legend.title = element_text(face = "bold", size = 14))+theme(axis.title.x = element_text(face="bold",size = 14, color="black"), axis.text.x = element_text(face="bold",size = 14, color="black",angle = 45, vjust = 0.5, hjust = 0.5), axis.title.y = element_text(face="bold",size = 14, color="black"), axis.text.y = element_text(face="bold", size = 14, color="black"))






#####################################################################
#####################################################################
#####################################################################
#Legends

p_tmp <- p1+  theme(strip.text.x = element_text(face = "bold", size = 18),strip.text.y = element_text(face = "bold", size = 18))+theme(legend.position="bottom")+  theme(legend.text = element_text(face = "bold", size = 18),legend.title = element_text(face = "bold", size = 18))+theme(axis.title.x = element_text(face="bold",size = 18, color="black"), axis.text.x = element_text(face="bold",size = 18, color="black",angle = 45, vjust = 0.5, hjust = 0.5), axis.title.y = element_text(face="bold",size = 18, color="black"), axis.text.y = element_text(face="bold", size = 18, color="black"))

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

pg3$heights = pg4$heights

p <- ggarrange(ggarrange(p1, p2, 
                         ncol = 2, labels = c("a", "b"),
                         widths = c(0.5,0.5)) , 
               ggarrange(pg3, pg4,
                         ncol = 2, labels = c("c", "d"),
                         widths = c(0.5,0.5)) ,
               nrow = 2, 
               labels = c(NA, NA),
               heights = c(0.52,0.48))

p <- cowplot::plot_grid( legend_b, p,  ncol = 1, rel_heights = c(.05, 1))

ggsave(filename=paste0("Figure3_UKB_09262024.pdf"),
       plot=p, 
       path="./family/manuscript/plots/final/",
       width=360, height=350, units="mm", dpi=320,device = cairo_pdf)

ggsave(filename=paste0("Figure3_UKB_09262024.png"),
       plot=p, 
       path="./family/manuscript/plots/final/",
       width=360, height=350, units="mm", dpi=320)

