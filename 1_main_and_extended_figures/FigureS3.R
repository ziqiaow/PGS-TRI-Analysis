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

#plot bias and sd
p <- ggplot(data_plot[which(data_plot$true_beta==0 & data_plot$true_delta != -0.2),], aes(x = Family_Size, bias,
                                                                                     ymin = bias - sd,
                                                                                     ymax = bias + sd,colour = Method , fill = Method)) +
  scale_x_discrete('# of Families') + scale_y_continuous('Bias +- SD \n',limits = c(-0.4,0.5))+
  geom_errorbar(position = position_dodge(width=0.6),size=1,width=0.4, aes(color=Method))+ scale_color_npg()  + 
  geom_hline(yintercept = 0, colour = "grey", linetype='dashed', size=0.5) +
  geom_point(shape = 19, size = 2,position = position_dodge(width=0.6)) +
  scale_color_manual(values = pal_nejm("default")(5)[c(1,2)]) +
  facet_grid(cor~truedelta)+ theme_bw()+      
  theme(plot.title = element_text(color = "black", size = 7, face = "bold")) +                                                          
  My_Theme+theme(legend.position="top",legend.text = element_text(face = "bold", size = 7),legend.title = element_text(face = "bold", size = 7)) 

ggsave(filename=paste0("Extended_Figure3.pdf"),
       plot=p, 
       path="./family/manuscript/draft4/plots/extended/",
       width=180, height=160, units="mm", dpi=400,device = cairo_pdf)
ggsave(filename=paste0("Extended_Figure3.png"),
       plot=p, device="png",
       path="./family/manuscript/draft4/plots/extended/",
       width=180, height=160, units="mm", dpi=400)



