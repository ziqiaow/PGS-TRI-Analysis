#----------------------------------------------------------
#Extended Figure 3
#----------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggsci)
load("./family/simulation/nurture/redo/08222024/results_summary_pgsmain_misfit.rda")
data_plot=bind_rows(summary_beta)
data_plot$cor=c(rep("Cor(\u03b1,\u03bc) = 0.25",24*2),rep("Cor(\u03b1,\u03bc) = 0",24))
data_plot$truedelta=rep(c(rep("\u03B4-IDE = 0",8),rep("\u03B4-IDE = 0.2",8),rep("\u03B4-IDE = -0.2",8)),3)
data_plot$truebeta=c(rep("\u03B2(PGS) = 0",24),rep("\u03B2(PGS) = 0.4",24),rep("\u03B2(PGS) = 0",24))
data_plot$true_beta=c(rep(0,24),rep(0.4,24),rep(0,24))
data_plot$true_delta=rep(c(rep(0,8),rep(0.2,8),rep(-0.2,8)),3)

data_plot$cor=factor(data_plot$cor,levels=c("Cor(\u03b1,\u03bc) = 0","Cor(\u03b1,\u03bc) = 0.25"))
#data_plot$cor=factor(data_plot$cor,levels=c("No population stratification","Cor(\u03b1,\u03bc) = 0.25"))
data_plot$truedelta=factor(data_plot$truedelta,levels=unique(data_plot$truedelta))
colnames(data_plot)[7]="Method"
colnames(data_plot)[8]="Family_Size"
data_plot$Family_Size=as.character(data_plot$Family_Size)
data_plot$Family_Size=factor(data_plot$Family_Size,level=c("200","500","1000","2000"))
data_plot$Method[which(data_plot$Method=="fam")]="PGS-TRI"
data_plot$Method[which(data_plot$Method=="logistic")]="Logistic Regression"
data_plot$Method=factor(data_plot$Method,level=c("PGS-TRI","Logistic Regression"))




#plot bias and sd
p <- ggplot(data_plot[which(data_plot$true_beta==0 & data_plot$true_delta != -0.2),], aes(x = Family_Size, bias,
                                                                                     ymin = bias - sd,
                                                                                     ymax = bias + sd,colour = Method , fill = Method)) +
  scale_x_discrete('# of Families') + scale_y_continuous('Bias +- SD \n',limits = c(-0.4,0.5))+
  geom_errorbar(position = position_dodge(width=0.6),size=2,width=0.4, aes(color=Method))+ scale_color_npg()  + 
  geom_hline(yintercept=0, col = "grey")  + 
  geom_point(shape = 19, size = 5,position = position_dodge(width=0.6)) +
  # theme(legend.position="none") + 
  scale_color_manual(values = pal_nejm("default")(5)[c(1,2)]) +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12))+ #+ labs(title="Simulation Investigation of the PGS Main Effect in Case-Parent Trio Study") +
  theme(plot.title = element_text(color = "black", size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  facet_grid(cor~truedelta)+ theme_bw()+                                                                # Change font size
  theme(strip.text.x = element_text(face = "bold", size = 18),strip.text.y = element_text(face = "bold", size = 18))+theme(legend.position="top")+  theme(legend.text = element_text(face = "bold", size=18),legend.title = element_text(face = "bold", size=18))+theme(axis.title.x = element_text(face="bold",size=18, color="black"), axis.text.x = element_text(face="bold",size=18, color="black"), axis.title.y = element_text(face="bold",size=18, color="black"), axis.text.y = element_text(face="bold", size=18, color="black"))

ggsave(filename=paste0("Extended_Figure3_09262024.pdf"),
       plot=p, 
       path="./family/manuscript/plots/final/",
       width=12,height=10, dpi=320,device = cairo_pdf)
ggsave(filename=paste0("Extended_Figure3_09262024.png"),
       plot=p, device="png",
       path="./family/manuscript/plots/final/",
       width=12,height=10, dpi=320)



