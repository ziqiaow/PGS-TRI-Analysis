

#----------------------------------------------------------
#Figure 3
#----------------------------------------------------------
setwd("./family/manuscript/draft4/plots/main/")
path<-"./family/manuscript/draft4/plots/main/"

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
  text = element_text(color = "black",size = 6),
  axis.text = element_text(color = "black", size = 7),
  axis.title = element_text(size = 7),
  axis.text.x = element_text(angle = 45,  vjust = 0.55, hjust = 0.5),
  strip.text.x = element_text(size = 6),strip.text.y = element_text(size = 6),
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

load("./family/simulation/UKB/match_geo/new_score/clump/PS/PGS_main/results_summary_multichild.rda")

data_plot=bind_rows(summary_type1error_power[3],summary_type1error_power[4])
data_plot = data_plot[-which(data_plot$method == "pcfather" | data_plot$method == "pcmother"| data_plot$method == "logistic_pcfamo" | data_plot$method == "logistic_pc10"),]
data_plot$Family_Size = c(rep("#Trios = 1000",5),rep("#Trios = 2000",5))

data_plot$method = rep(c("PGS-TRI","Logistic Regression","pTDT","Logistic Regression + PC10 Geo","Logistic Regression + PGS(MF)"),2)
data_plot$method=factor(data_plot$method,levels=c("PGS-TRI","Logistic Regression","Logistic Regression + PC10 Geo","Logistic Regression + PGS(MF)","pTDT"))

colnames(data_plot)[2]="Method"
data_plot$measure = "PS"
data_plot1 = data_plot

#load data with population stratification bias and assortative mating matched by geographical regions and educational attainments
load("./family/simulation/UKB/new_score/clump/PS_AM/PGS_main/results_summary_multichild.rda")
#type I error
data_plot=bind_rows(summary_type1error_power[3],summary_type1error_power[4])
data_plot = data_plot[-which(data_plot$method == "pcfather" | data_plot$method == "pcmother"| data_plot$method == "logistic_pcfamo" | data_plot$method == "logistic_pc10"),]
data_plot$Family_Size = c(rep("#Trios = 1000",5),rep("#Trios = 2000",5))

data_plot$method = rep(c("PGS-TRI","Logistic Regression","pTDT","Logistic Regression + PC10 Geo","Logistic Regression + PGS(MF)"),2)
data_plot$method=factor(data_plot$method,levels=c("PGS-TRI","Logistic Regression","Logistic Regression + PC10 Geo","Logistic Regression + PGS(MF)","pTDT"))


colnames(data_plot)[2]="Method"
data_plot$measure = "PS + AM"


data_plot = rbind(data_plot1,data_plot)
data_plot = data_plot[which(data_plot$Family_Size == "#Trios = 1000"),]
data_plot$measure = factor(data_plot$measure,levels = unique(data_plot$measure))
data_plot$effect = "DE"

p1 <- ggplot(data=data_plot, aes(x=Method, y=type1error_beta_prs , fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width = 0.9,linewidth = 0.2)+
  scale_fill_manual(values = c(pal_nejm("default")(5)[c(1)],"Royal Blue","Cornflower Blue","Light Sky Blue",pal_nejm("default")(5)[c(4)])) +
  scale_x_discrete(NULL) +
  scale_y_continuous('Type I Error',breaks = c(0.05,0.5,1),limits = c(0,1))+
  geom_hline(yintercept=0.05, col = "coral",linetype = "dashed",size=0.5)  + 
  facet_grid(measure~effect)+
  theme_bw()+      
  My_Theme+theme(legend.position="") 


#####################################################################
#####################################################################
#####################################################################
# Panel C


#plot bias and sd for UKB simulation results
#load data with population stratification bias matched by geographical regions only

load("./family/simulation/UKB/match_geo/new_score/clump/PS/PGS_main/results_summary_multichild.rda")


#plot bias and sd

data_plot=bind_rows(summary_bias[c(1,2,3,4)]) 
data_plot = data_plot[-which(data_plot$method == "pcfather" | data_plot$method == "pcmother"| data_plot$method == "logistic_pcfamo" | data_plot$method == "logistic_pc10"),]
data_plot$truebeta=c(rep("\u03B2(G) = 0.4",10),rep("\u03B2(G) = 0",10))
data_plot$Family_Size = c(rep("#Trios = 1000",5),rep("#Trios = 2000",5),rep("#Trios = 1000",5),rep("#Trios = 2000",5))
data_plot$truebeta=factor(data_plot$truebeta,levels=c("\u03B2(G) = 0","\u03B2(G) = 0.4"))
data_plot$method = rep(c("PGS-TRI","Logistic Regression","pTDT","Logistic Regression + PC10 Geo","Logistic Regression + PGS(MF)"),4)
data_plot$method=factor(data_plot$method,levels=c("PGS-TRI","Logistic Regression","Logistic Regression + PC10 Geo","Logistic Regression + PGS(MF)","pTDT"))


colnames(data_plot)[2]="Method"

data_plot3=bind_rows(summary_empirical_sd[c(1,2,3,4)])
data_plot3 = data_plot3[-which(data_plot3$method == "pcfather" | data_plot3$method == "pcmother"| data_plot3$method == "logistic_pcfamo" | data_plot3$method == "logistic_pc10"),]

data_plot=cbind(data_plot,data_plot3[,1])
colnames(data_plot)[5]=c("sd1")
data_plot$measure = "PS"

data_plot1 = data_plot


#load data with population stratification bias and assortative mating matched by geographical regions and educational attainments
load("./family/simulation/UKB/new_score/clump/PS_AM/PGS_main/results_summary_multichild.rda")

#plot bias and sd

data_plot=bind_rows(summary_bias[c(1,2,3,4)])
data_plot = data_plot[-which(data_plot$method == "pcfather" | data_plot$method == "pcmother"| data_plot$method == "logistic_pcfamo" | data_plot$method == "logistic_pc10"),]
data_plot$truebeta=c(rep("\u03B2(G) = 0.4",10),rep("\u03B2(G) = 0",10))
data_plot$Family_Size = c(rep("#Trios = 1000",5),rep("#Trios = 2000",5),rep("#Trios = 1000",5),rep("#Trios = 2000",5))
data_plot$truebeta=factor(data_plot$truebeta,levels=c("\u03B2(G) = 0","\u03B2(G) = 0.4"))
data_plot$method = rep(c("PGS-TRI","Logistic Regression","pTDT","Logistic Regression + PC10 Geo","Logistic Regression + PGS(MF)"),4)
data_plot$method=factor(data_plot$method,levels=c("PGS-TRI","Logistic Regression","Logistic Regression + PC10 Geo","Logistic Regression + PGS(MF)","pTDT"))

colnames(data_plot)[2]="Method"

data_plot3=bind_rows(summary_empirical_sd[c(1,2,3,4)]) 
data_plot3 = data_plot3[-which(data_plot3$method == "pcfather" | data_plot3$method == "pcmother"| data_plot3$method == "logistic_pcfamo" | data_plot3$method == "logistic_pc10"),]

data_plot=cbind(data_plot,data_plot3[,1])
colnames(data_plot)[5]=c("sd1")
data_plot$measure = "PS + AM"


data_plot = rbind(data_plot1,data_plot)
data_plot = data_plot[-which(data_plot$Method == "pTDT"),]
data_plot$Family_Size = factor(data_plot$Family_Size)


data_plot = data_plot[which(data_plot$Family_Size=="#Trios = 1000"),]
data_plot$measure = factor(data_plot$measure,levels = unique(data_plot$measure))
data_plot$effect = "DE"
p3 <- ggplot(data_plot, aes(x = truebeta, bias_beta_prs,
                            ymin = bias_beta_prs - sd1,
                            ymax = bias_beta_prs + sd1,colour = Method , fill = Method)) +
  scale_x_discrete(NULL) +
  scale_y_continuous('Bias +- SD \n',lim = c(-0.6, 0.6)
  ) +
  geom_errorbar(position = position_dodge(width=0.6),size=0.5,width=0.4, aes(color=Method))+ 
  geom_hline(yintercept = 0, colour = "grey", linetype='dashed', size=0.5) +
  scale_color_manual(values = c(pal_nejm("default")(5)[c(1)],"Royal Blue","Cornflower Blue","Light Sky Blue",pal_nejm("default")(5)[c(4)])) +
  geom_point(shape = 19, size = 1,position = position_dodge(width=0.6))+
  facet_grid(measure~effect)+ 
  theme_bw()+      
  My_Theme+theme(legend.position="") 


#####################################################################
#####################################################################
#####################################################################
# Panel B

#type I error for UKB simulation results for IDE
#load data with PS matched by geographical regions only
load("./family/simulation/UKB/match_geo/new_score/clump/PS/nurture/multichild/results_summary.rda")

#type I error
data_plot=bind_rows(summary_delta)
data_plot = data_plot[-which(data_plot$method == "logistic_pc10" | data_plot$method == "logistic_pc10_geo"),]

data_plot$truedelta=rep(c(rep("\u03B4-IDE = 0",3),rep("\u03B4-IDE = 0.2",3),rep("\u03B4-IDE = -0.2",3)),2)
data_plot$type = "\u03B4-IDE"
data_plot$truedelta=factor(data_plot$truedelta,levels=unique(data_plot$truedelta))
colnames(data_plot)[7]="Method"
colnames(data_plot)[9]="Family_Size"
data_plot$Family_Size=as.character(data_plot$Family_Size)
data_plot$Family_Size=factor(data_plot$Family_Size,level=c("1000","2000"))
data_plot$Method = c(rep(c("PGS-TRI","Logistic Regression + PGS(MF)","PGS-TRI-centered"),6))
data_plot$Method=factor(data_plot$Method,levels=c("PGS-TRI","PGS-TRI-centered","Logistic Regression + PGS(MF)"))
data_plot$coverage_prob=data_plot$coverage_prob*100

#Type I error
data_plot1=data_plot[which(data_plot$true_delta==0 ),]
data_plot1$measure = "PS"


#load data with population stratification bias and assortative mating matched by geographical regions and educational attainments
load("./family/simulation/UKB/new_score/clump/PS_AM/nurture/multichild/results_summary.rda")

#type I error
data_plot=bind_rows(summary_delta)
data_plot = data_plot[-which(data_plot$method == "logistic_pc10" | data_plot$method == "logistic_pc10_geo"),]
data_plot$truedelta=rep(c(rep("\u03B4-IDE = 0",3),rep("\u03B4-IDE = 0.2",3),rep("\u03B4-IDE = -0.2",3)),2)
data_plot$type = "\u03B4-IDE"

data_plot$truedelta=factor(data_plot$truedelta,levels=unique(data_plot$truedelta))
colnames(data_plot)[7]="Method"
colnames(data_plot)[9]="Family_Size"
data_plot$Family_Size=as.character(data_plot$Family_Size)
data_plot$Family_Size=factor(data_plot$Family_Size,level=c("1000","2000"))
data_plot$Method = c(rep(c("PGS-TRI","Logistic Regression + PGS(MF)","PGS-TRI-centered"),6))
data_plot$Method=factor(data_plot$Method,levels=c("PGS-TRI","PGS-TRI-centered","Logistic Regression + PGS(MF)"))
data_plot$coverage_prob=data_plot$coverage_prob*100

#Type I error
data_plot2=data_plot[which(data_plot$true_delta==0 ),]
data_plot2$measure = "PS + AM"


data_plot = rbind(data_plot1,data_plot2)
data_plot$measure = factor(data_plot$measure,levels = unique(data_plot$measure))

data_plot = data_plot[which(data_plot$Family_Size == "1000"),]


p2 <- ggplot(data=data_plot, aes(x=Method, y=type1error.power , fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width = 0.9,linewidth=0.2)+
  scale_fill_manual(values = c(pal_nejm("default")(5)[c(1,3)],"Light Sky Blue")) +
  scale_x_discrete(NULL) +
  scale_y_continuous('Type I Error',breaks = c(0.05,0.5),limits = c(0,0.5))+
  geom_hline(yintercept=0.05, col = "coral",linetype = "dashed",size=0.5)  + 
  facet_grid(measure~type,scales = "free")+ 
  theme_bw()+      
  My_Theme+theme(legend.position="")

#####################################################################
#####################################################################
#####################################################################
# Panel D
#load data with PS matched by geographical regions only
load("./family/simulation/UKB/match_geo/new_score/clump/PS/nurture/multichild/results_summary.rda")
data_plot=bind_rows(summary_delta)
data_plot = data_plot[-which(data_plot$method == "logistic_pc10" | data_plot$method == "logistic_pc10_geo"),]
data_plot$truedelta=rep(c(rep("\u03B4-IDE = 0",3),rep("\u03B4-IDE = 0.2",3),rep("\u03B4-IDE = -0.2",3)),2)

data_plot$type = "\u03B4-IDE"
data_plot$truedelta=factor(data_plot$truedelta,levels=unique(data_plot$truedelta))
colnames(data_plot)[7]="Method"
colnames(data_plot)[9]="Family_Size"
data_plot$Family_Size=as.character(data_plot$Family_Size)
data_plot$Family_Size=factor(data_plot$Family_Size,level=c("1000","2000"))
data_plot$Method = c(rep(c("PGS-TRI","Logistic Regression + PGS(MF)","PGS-TRI-centered"),6))
data_plot$Method=factor(data_plot$Method,levels=c("PGS-TRI","PGS-TRI-centered","Logistic Regression + PGS(MF)"))
data_plot$coverage_prob=data_plot$coverage_prob*100
data_plot1 = data_plot
data_plot1$measure = "PS"


#load data with population stratification bias and assortative mating matched by geographical regions and educational attainments
load("./family/simulation/UKB/new_score/clump/PS_AM/nurture/multichild/results_summary.rda")
data_plot=bind_rows(summary_delta)
data_plot = data_plot[-which(data_plot$method == "logistic_pc10" | data_plot$method == "logistic_pc10_geo"),]
data_plot$truedelta=rep(c(rep("\u03B4-IDE = 0",3),rep("\u03B4-IDE = 0.2",3),rep("\u03B4-IDE = -0.2",3)),2)
data_plot$type = "\u03B4-IDE"
data_plot$truedelta=factor(data_plot$truedelta,levels=unique(data_plot$truedelta))
colnames(data_plot)[7]="Method"
colnames(data_plot)[9]="Family_Size"
data_plot$Family_Size=as.character(data_plot$Family_Size)
data_plot$Family_Size=factor(data_plot$Family_Size,level=c("1000","2000"))
data_plot$Method = c(rep(c("PGS-TRI","Logistic Regression + PGS(MF)","PGS-TRI-centered"),6))
data_plot$Method=factor(data_plot$Method,levels=c("PGS-TRI","PGS-TRI-centered","Logistic Regression + PGS(MF)"))
data_plot$coverage_prob=data_plot$coverage_prob*100
data_plot2 = data_plot
data_plot2$measure = "PS + AM"




data_plot = rbind(data_plot1,data_plot2)
data_plot$measure = factor(data_plot$measure,levels = unique(data_plot$measure))

data_plot = data_plot[which(data_plot$Family_Size == "1000"),]
#plot bias and sd
p4 <- ggplot(data_plot, aes(x = truedelta, bias,
                      ymin = bias - sd,
                      ymax = bias + sd,colour = Method , fill = Method)) +
  scale_x_discrete('') +
  scale_y_continuous('Bias +- SD \n',lim = c(-0.6, 0.6)
  ) +
 geom_errorbar(position = position_dodge(width=0.6),size=0.5,width=0.4, aes(color=Method))+ 
  geom_hline(yintercept = 0, colour = "grey", linetype='dashed', size=0.5) +
  geom_point(shape = 19, size = 1,position = position_dodge(width=0.6)) +
  scale_color_manual(values =  c(pal_nejm("default")(5)[c(1,3)],"Light Sky Blue",pal_nejm("default")(5)[c(4)]))+
  facet_grid(measure~type,scales = "free")+ 
  theme_bw()+      
  My_Theme+theme(legend.position="")




#####################################################################
#####################################################################
#####################################################################
#Legends
data_plot$Method = as.character(data_plot$Method)
data_plot$Method[4] = "Logistic Regression + PC10 Geo"
data_plot$Method[5] = "Logistic Regression"
data_plot$Method[6] = "pTDT"

data_plot$Method=factor(data_plot$Method,levels=c("PGS-TRI","PGS-TRI-centered","Logistic Regression","Logistic Regression + PC10 Geo","Logistic Regression + PGS(MF)","pTDT"))


p_tmp <- ggplot(data=data_plot, aes(x=truedelta, y=type1error.power , fill=Method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width = 0.2,linewidth = 0.2)+
  scale_fill_manual(values = c(pal_nejm("default")(5)[c(1,3)],"Royal Blue","Cornflower Blue","Light Sky Blue",pal_nejm("default")(5)[c(4)])) +
  scale_x_discrete(NULL) +
  scale_y_continuous(NULL,breaks = c(0.05,0.5,1),limits = c(0,0.5))

p_tmp <- p_tmp+ theme(legend.text = element_text( size = 6),legend.title = element_text( size = 6),legend.key.size = unit(3, 'mm'))+
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) 

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
                         font.label = list(size = 7),
                         widths = c(0.5,0.5)) , 
               ggarrange(pg3, pg4,
                         ncol = 2, labels = c("c", "d"),
                         font.label = list(size = 7),
                         widths = c(0.5,0.5)) ,
               nrow = 2, 
               labels = c(NA, NA),
               heights = c(0.52,0.48))

p <- cowplot::plot_grid( legend_b, p,  ncol = 1, rel_heights = c(.05, 1))

ggsave(filename=paste0("Figure3.pdf"),
       plot=p, 
       path="./family/manuscript/draft4/plots/main/",
       width=180, height=180, units="mm", dpi=1000,device = cairo_pdf)

ggsave(filename=paste0("Figure3.png"),
       plot=p, 
       path="./family/manuscript/draft4/plots/main/",
       width=180, height=180, units="mm", dpi=1000)

