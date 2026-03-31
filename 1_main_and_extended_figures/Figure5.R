
#----------------------------------------------------------
#Figure 5
#----------------------------------------------------------

library(patchwork)
library(dplyr, warn = FALSE)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

# Panel A
source("./family/R/PGS.TRI/R/PGS-TRI.R")
source("./family/R/PGS.TRI/R/pTDT.R")

#load data and data cleaning
library(data.table)
library(fastDummies)
library(ggplot2)
library(ggsci)
library(dplyr)

load("./family/data/autism/v3/PRS_SPARK_PCstandardized.RData")
load("./family/data/autism/v3/roles_id_complete_casetrio_ancestry.RData")

#check reference data PC mean
load("./family/data/1000G_HGDP/PCA/pc_1KG_race.RData")

plot_c = dat[match(id_final$subject_sp_id,dat$IID),c(1,7,11:20)]
plot_c = plot_c[which(plot_c$superclass %in% c("AFR", "AMR", "EUR", "EAS", "SAS")),]

#plot together with 1000G+HGDP 
#combine both dataset together
pc_1KG_race$Cohort = "1000G+HGDP"
colnames(pc_1KG_race)[16] = "superclass"
plot_c$Cohort = "SPARK"

pc_1KG_spark = rbind(pc_1KG_race[,c("IID","PC1","PC2","superclass","Cohort")],plot_c[,c("IID","PC1","PC2","superclass","Cohort")])


#load Finnish people
fin_ID = fread("./family/data/1000G_HGDP/PCA/sample_fin.tsv")
id = intersect(fin_ID$`Sample name`,pc_1KG_race$IID)
pc_1KG_race = pc_1KG_race[match(id,pc_1KG_race$IID),]
pc_center = apply(pc_1KG_race[,c(5:14)],2,mean)

#combine both dataset together
colnames(pc_1KG_race)[16] = "superclass"
pc_1KG_race$superclass = "Finnish EUR (1000G+HGDP)"
pc_1KG_spark = rbind(pc_1KG_race[,c("IID","PC1","PC2","PC3","superclass","Cohort")],plot_c[,c("IID","PC1","PC2","PC3","superclass","Cohort")])
pc_1KG_spark$superclass = factor(pc_1KG_spark$superclass,levels = c("EUR","AFR","AMR","EAS","SAS","Finnish EUR (1000G+HGDP)"))


col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

p1 <- ggplot(pc_1KG_spark, aes(x = PC1, y = PC2))+#, color = factor(superclass), shape = factor(Cohort), size = factor(Cohort))) + 
  geom_point(data = subset(pc_1KG_spark, Cohort != "1000G+HGDP"),
             aes(color = superclass), alpha = 1, size = 0.5) +  
  #scale_color_jco() +
  geom_point(data = subset(pc_1KG_spark, Cohort == "1000G+HGDP"),
             aes(color = superclass),  alpha = 1, size = 0.5) +  
  scale_color_manual(values=c("#DE3163",pal_jco("default")(10)[c(1,2,3,5)],"dark red"),labels = c("EUR","AFR","AMR","EAS","SAS","Finnish EUR (1000G+HGDP)")) +
  ggforce::geom_circle(
    data = subset(pc_1KG_spark, Cohort == "1000G+HGDP") %>%   
      summarise(
        x0 = mean(PC1),  
        y0 = mean(PC2),  
        r = 0.01         
      ),
    aes(x0 = x0, y0 = y0, r = r),
    color = "black",
    # linetype = "dashed",
    linewidth = 0.8,
    inherit.aes = FALSE  
  ) +
  coord_fixed()+ theme_bw()+ theme(
    panel.background = element_blank(), strip.background = element_rect(colour = NA, fill = NA),
    strip.text.y = element_text(face = "bold",size = 7, color = "black"),
    strip.text = element_text(face = "bold",size = 7, color = "black"),
    panel.border = element_rect(fill = NA, color = "black",size=0.15),
    legend.position = "bottom",
    axis.text = element_text(face = "bold",size = 7, color = "black"),
    axis.title = element_text(face = "bold",size = 7, color = "black",hjust = 0.5),
    plot.title = element_text(face = "bold", hjust = 0.5,size = 7, color = "black"),
    legend.title=element_text(face = "bold", size=7), 
    legend.text=element_text(face = "bold", size=7),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.15),
    axis.ticks = element_line(size = 0.15)
  )+
  theme(
    plot.margin = margin(0, 1, 0, 1, "mm")
  ) + theme(legend.key.size = unit(4, 'mm'))+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE),color = guide_legend(override.aes = list(size = 3),title="")) 



p2 <- ggplot(pc_1KG_spark, aes(x = PC1, y = PC3))+#, color = factor(superclass), shape = factor(Cohort), size = factor(Cohort))) + 
  geom_point(data = subset(pc_1KG_spark, Cohort != "1000G+HGDP"),
             aes(color = superclass), alpha = 1, size = 0.5) +  
  #scale_color_jco() +
  geom_point(data = subset(pc_1KG_spark, Cohort == "1000G+HGDP"),
             aes(color = superclass),  alpha = 1, size = 0.5) +  
  scale_color_manual(values=c("#DE3163",pal_jco("default")(10)[c(1,2,3,5)],"dark red"),labels = c("EUR","AFR","AMR","EAS","SAS","Finnish EUR (1000G+HGDP)")) +
  ggforce::geom_circle(
    data = subset(pc_1KG_spark, Cohort == "1000G+HGDP") %>%   
      summarise(
        x0 = mean(PC1),  
        y0 = mean(PC3),  
        r = 0.01         
      ),
    aes(x0 = x0, y0 = y0, r = r),
    color = "black",
    # linetype = "dashed",
    linewidth = 0.8,
    inherit.aes = FALSE  
  ) +
  coord_fixed()+ theme_bw()+guides(color = guide_legend(title=""))+ theme(
    panel.background = element_blank(), strip.background = element_rect(colour = NA, fill = NA),
    strip.text.y = element_text(face = "bold",size = 6, color = "black"),
    strip.text = element_text(face = "bold",size = 6, color = "black"),
    panel.border = element_rect(fill = NA, color = "black",size=0.15),
    legend.position = "bottom",
    axis.text = element_text(face = "bold",size = 7, color = "black"),
    axis.title = element_text(face = "bold",size = 7, color = "black",hjust = 0.5),
    plot.title = element_text(face = "bold", hjust = 0.5,size = 7, color = "black"),
    legend.title=element_text(face = "bold", size=7), 
    legend.text=element_text(face = "bold", size=7),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.15),
    axis.ticks = element_line(size = 0.15)
  )+
  theme(
    plot.margin = margin(0, 1, 0, 1, "mm")
  )+ theme(legend.key.size = unit(4, 'mm'))+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE),color = guide_legend(override.aes = list(size = 3),title="")) 



# Panel B
#Make a figure for genetic distance and RR
#load("./family/data/autism/v3/GxE/PGS_GD_plot.RData")
load("./family/data/autism/v3/GxE/PGS_GD_plot_PCmeanstandardizeonly_raw.RData")

#main effect
# Reduce the opacity of the grid lines: Default is 255
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

p3 <- ggplot(
  data = table_all[[2]],
  aes(x = GD, y = OR, ymin = ci_low, ymax = ci_high)
)+#ylim(0.8,1.35) +
  geom_pointrange(aes(col = GD_group_ancestry,shape = sig), position=position_dodge(width=0.8),size = 0.3) +scale_shape_manual(values = c(8, 16),name="P < 0.05")+
  geom_hline(yintercept = 1, colour = "grey", linetype='dashed', size=0.5) +
  xlab("Mean Genetic Distance from Reference Dataset") +
  ylab("Relative Risk for ASD (95% CI)") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, col = GD_group_ancestry), width = 0.2, cex = 0.4, position=position_dodge(width=0.5)) +
  # scale_color_manual(values = rep(c("steel blue","dark red"), 1),name="Genetic Effect") + #"orange","purple"
  # facet_wrap(~race, strip.position = "top", nrow = 2, scales = "free_y") +
  scale_color_ucscgb() + theme_classic() +
  scale_y_continuous(breaks = c(0.8, 1, 1.2,1.4))+
  theme(
    panel.background = element_blank(), strip.background = element_rect(colour = NA, fill = NA),
    strip.text.y = element_text(face="bold",size = 7, color = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.text = element_text(face="bold",size = 7, color = "black"),
    panel.border = element_rect(fill = NA, color = "black",size=0.15),
    legend.position = "right",
    axis.text = element_text(face="bold",size = 7, color = "black"),
    axis.title = element_text(face="bold",size = 7, color = "black",hjust = 0.5),
    legend.title=element_text(face="bold",size=7), 
    legend.text=element_text(size=5),
    legend.spacing.y = unit(0, "mm"),
    legend.key.spacing.y = unit(0, "mm"),
    legend.key.size = unit(4,"mm"),
   # panel.border = element_rect(size = 0.15),
   # strip.background = element_rect(size = 0.15),
    axis.line = element_line(size = 0.15),
    axis.ticks = element_line(size = 0.15)
    
  ) +
  annotate(
    "text",
    x = -Inf,          
    y = -Inf,          
    label = "PGS x Genetic Distance Interaction:\nlog RR (SD) = -0.04 (0.01) P < 1e-04",  
    fontface =2,
    hjust = -0.1,      
    vjust = -0.5,      
    color = "black",    
    size = 2     
  )+
  theme(
    plot.margin = margin(4, 4, 4, 4, "mm")
  )# + guides(color = guide_legend(reverse = TRUE))  


###############################################################
###############################################################
###############################################################


p <- ggarrange(ggarrange(p1,p2,common.legend = T,legend="top" ,align = "hv",labels = c("a", "b"),font.label = list(size = 7),ncol = 2,widths = c(0.5,0.5)),p3,
               nrow = 2, 
               labels = c(NA, "c"),
               font.label = list(size = 7),
               heights = c(0.55,0.45)
)

ggsave(filename=paste0("Figure5.pdf"),
       plot=p, 
       path="./family/manuscript/draft4/plots/main/",
       width=180, height=170, units="mm", dpi=1000,device = cairo_pdf)

ggsave(filename=paste0("Figure5.png"),
       plot=p, device="png",
       path="./family/manuscript/draft4/plots/main/",
       width=180, height=170, units="mm", dpi=1000)



