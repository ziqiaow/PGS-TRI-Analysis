

#----------------------------------------------------------
#Figure 4
#----------------------------------------------------------
library(patchwork)
library(dplyr, warn = FALSE)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(data.table)
library(viridis)
library(ggsci)
library(RColorBrewer)

# Panel A
# Forest plot for main effect of Autism, and other traits

setwd("./family/data/autism/v3")
load("./family/data/autism/v3/results_ASD_all.RData")

main=data.frame(direct_effect_homo[which(direct_effect_homo$standardize == "PC based (with 1KG fixed SD)"),])
main = main[c(2:6,1),]
main$race=c("European","African","Americas",
            "East Asian","South Asian","Cross Population")

nurture=data.frame(indirect_effect_homo[which(indirect_effect_homo$standardize == "PC based (with 1KG fixed SD)"),])
nurture = nurture[c(2:6,1),]
nurture$race=c("European","African","Americas",
               "East Asian","South Asian","Cross Population")

table_asd=main

table_asd$effect = c("DE")



table_asd$OR=exp(table_asd$Estimate)
table_asd$ci_high=exp(table_asd$Estimate+qnorm(0.975)*table_asd$Std.Error)
table_asd$ci_low=exp(table_asd$Estimate-qnorm(0.975)*table_asd$Std.Error)


table_plot = table_asd
table_plot$race = factor(table_plot$race, levels = rev(unique(table_plot$race)))

table_plot1= table_plot

table_asd=nurture

table_asd$effect = c("\u03B4-IDE")




table_asd$OR=exp(table_asd$Estimate)
table_asd$ci_high=exp(table_asd$Estimate+qnorm(0.975)*table_asd$Std.Error)
table_asd$ci_low=exp(table_asd$Estimate-qnorm(0.975)*table_asd$Std.Error)


table_plot = table_asd
table_plot$race = factor(table_plot$race, levels = rev(unique(table_plot$race)))

table_plot_final = rbind(table_plot,table_plot1)
table_plot_final$sig = "no"
table_plot_final$sig[which(table_plot_final$Pvalue < 0.05)] = "yes"
table_plot_final$effect = factor(table_plot_final$effect,levels = rev(c("DE","\u03B4-IDE")))


#main effect
# Reduce the opacity of the grid lines: Default is 255
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)
table_plot_final$trait = "PGS for ASD"
p1<-ggplot(
  data = table_plot_final,
  aes(x = race, y = OR, ymin = ci_low, ymax = ci_high, group = effect)
)  +scale_y_continuous(limits = c(0.75,1.5),breaks = c(0.8, 1, 1.2,1.4))+
  geom_pointrange(aes(col = effect,shape = sig), position=position_dodge(width=0.5),
  fatten = 0) +  
  geom_point(aes(col = effect, shape = sig),
             position = position_dodge(width = 0.5),
             size = 1.2) +     
scale_shape_manual(values = c(16, 8))+
  geom_hline(yintercept = 1, colour = "grey", linetype='dashed', size=0.5) +
  xlab(NULL) +
  ylab("Relative Risk for ASD (95% CI)") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, col = effect),width = 0.4, cex = 0.5, position=position_dodge(width=0.5)) +
  scale_color_manual(values = rep(c("steel blue","dark red"), 1)) + #"orange","purple"
theme_classic() +
  facet_wrap(~trait, strip.position = "top", nrow = 2, scales = "free_y") +
  theme(
    panel.background = element_blank(), strip.background = element_rect(colour = NA, fill = NA),
   # strip.background = element_rect(size = 0.15),
    axis.line = element_line(size = 0.15),
    axis.ticks = element_line(size = 0.15),
    strip.text.y = element_text(size = 7, color = "black"),
    #panel.grid.major.y = element_line(colour = col_grid, size = 1 ),
   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7, color = "black"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none",
    legend.text=element_text(size=7),
    legend.title = element_text(size=7),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 7, color = "black",hjust = 0.5),
    plot.title = element_text( hjust = 0.5,size = 7, color = "black")
  ) + coord_flip()





# Panel B
# Forest plot for main effect of other traits for All population
#https://stackoverflow.com/questions/73250489/how-to-align-table-with-forest-plot-ggplot2

setwd("./family/data/autism/v3")
load("./family/data/pleiotropy/redo/results_pleiotropy_all.RData")

library(tidyverse)
direct_effect_homo_save <-direct_effect_homo %>% 
  bind_rows(.id = "trait")
indirect_effect_homo_save <-indirect_effect_homo %>% 
  bind_rows(.id = "trait")
indirect_effect_homo_save = indirect_effect_homo_save[which(indirect_effect_homo_save$standardize == "PC based (with 1KG fixed SD)"),]
direct_effect_homo_save = direct_effect_homo_save[which(direct_effect_homo_save$standardize == "PC based (with 1KG fixed SD)"),]
direct_effect_homo_save$ancestry = indirect_effect_homo_save$ancestry = rep(c("all","EUR","AFR","AMR","EAS","SAS"),11)
direct_effect_homo_save = direct_effect_homo_save[,c(1,2,3,5,7)]
indirect_effect_homo_save = indirect_effect_homo_save[,c(1,2,3,5,7)]

pattern = c("edu","schizophrenia","adhd","chronotype","insomnia","bipolar","bipolar1","bipolar2","depression","neuroticism","bmi_prive")
direct_effect_homo_save = direct_effect_homo_save[match(pattern,direct_effect_homo_save$trait),]
indirect_effect_homo_save = indirect_effect_homo_save[match(pattern,indirect_effect_homo_save$trait),]

res_plot=indirect_effect_homo_save
res_plot$race = "PGS for Other Traits (Cross Population)"
table_asd=res_plot[,c(2,3,4,6)]

table_asd$Trait = c("Educational Attainment","Schizophrenia","ADHD","Chronotype","Insomnia",
                    "Bipolar Disorder","Bipolar Disorder 1","Bipolar Disorder 2",
                    "Major Depression","Neuroticism","BMI")



table_asd$OR=exp(table_asd$Estimate)
table_asd$ci_high=exp(table_asd$Estimate+qnorm(0.975)*table_asd$Std.Error)
table_asd$ci_low=exp(table_asd$Estimate-qnorm(0.975)*table_asd$Std.Error)

table_plot = table_asd
table_plot$Trait = factor(table_plot$Trait, levels = rev(c("Educational Attainment","Schizophrenia","ADHD","Chronotype","Insomnia",
                                                           "Bipolar Disorder","Bipolar Disorder 1","Bipolar Disorder 2",
                                                           "Major Depression","Neuroticism","BMI")))
table_plot$sig = "No"
table_plot$sig[which(table_plot$Pvalue < 0.05)] = "Yes"

table_plot1 = table_plot



res_plot=direct_effect_homo_save
res_plot$race = "PGS for Other Traits (Cross Population)"
table_asd=res_plot[,c(2,3,4,6)]

table_asd$Trait = c("Educational Attainment","Schizophrenia","ADHD","Chronotype","Insomnia",
                    "Bipolar Disorder","Bipolar Disorder 1","Bipolar Disorder 2",
                    "Major Depression","Neuroticism","BMI")



table_asd$OR=exp(table_asd$Estimate)
table_asd$ci_high=exp(table_asd$Estimate+qnorm(0.975)*table_asd$Std.Error)
table_asd$ci_low=exp(table_asd$Estimate-qnorm(0.975)*table_asd$Std.Error)

table_plot = table_asd
table_plot$Trait = factor(table_plot$Trait, levels = rev(c("Educational Attainment","Schizophrenia","ADHD","Chronotype","Insomnia",
                                                           "Bipolar Disorder","Bipolar Disorder 1","Bipolar Disorder 2",
                                                           "Major Depression","Neuroticism","BMI")))

table_plot$sig = "No"
table_plot$sig[which(table_plot$Pvalue < 0.05)] = "Yes"

table_plot_total = rbind(table_plot,table_plot1)
table_plot_total$effect = c(rep("DE",11),rep("\u03B4-IDE",11))
table_plot_total$effect = factor(table_plot_total$effect,levels = rev(c("DE","\u03B4-IDE")))
#main effect
# Reduce the opacity of the grid lines: Default is 255
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

p2 <- ggplot(
  data = table_plot_total,
  aes(x = Trait, y = OR, ymin = ci_low, ymax = ci_high, group = effect)
)+
  scale_y_continuous(limits = c(0.75,1.5),breaks = c(0.8, 1, 1.2,1.4))+
  geom_pointrange(aes(col = effect,shape = sig), position=position_dodge(width=0.5),
                  fatten = 0) +  
  geom_point(aes(col = effect, shape = sig),
             position = position_dodge(width = 0.5),
             size = 1.2) +     
  scale_shape_manual(values = c(16, 8),name="P < 0.05")+
  geom_hline(yintercept = 1, colour = "grey", linetype='dashed', size=0.5) +
  xlab(NULL) +
  ylab("Relative Risk for ASD (95% CI)") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, col = effect), width = 0.4, cex = 0.5, position=position_dodge(width=0.5)) +
  scale_color_manual(values = rep(c("steel blue","dark red"), 1),name="Genetic Effect") + #"orange","purple"
  facet_wrap(~race, strip.position = "top", nrow = 2, scales = "free_y") +
  theme_classic() +
  theme(
    panel.background = element_blank(), strip.background = element_rect(colour = NA, fill = NA),
  #  strip.background = element_rect(size = 0.15),
    axis.line = element_line(size = 0.15),
    axis.ticks = element_line(size = 0.15),
    strip.text.y = element_text(size = 7, color = "black"),
    #panel.grid.major.y = element_line(colour = col_grid, size = 1 ),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7, color = "black"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "right",
    legend.text=element_text(size=7),
    legend.title = element_text(size=7),
    axis.text = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 7, color = "black",hjust = 0.5),
    plot.title = element_text( hjust = 0.5,size = 7, color = "black")
  ) + coord_flip() + 
  guides(color = guide_legend(reverse = TRUE),shape = guide_legend(reverse = TRUE))


# Panel B
# Miami plot for Transcriptomics PGS main results in Asian and EUR



library(readr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(latex2exp)

thres <- 0.05

My_Theme = theme(
  panel.background = element_blank(), 
  title = element_text(size = 7),
  text = element_text(size = 6)
)

###############################################################
###############################################################
###############################################################

# Panel C, PGS main
#load annotation

type="RNAseq"
direct_effect = fread(paste0("./family/data/phewas/results/homo/",type,"_direct_anno_gene.csv"))
indirect_effect = fread(paste0("./family/data/phewas/results/homo/",type,"_indirect_anno_gene.csv"))
identical(rownames(direct_effect),rownames(indirect_effect))

plot_man=direct_effect
plot_man=plot_man[,c(1,2,6:8,22:39)]
colnames(plot_man)[c(3:5)] = c("beta_direct","se_direct","p_direct")
plot_man$beta_indirect = indirect_effect$beta_EUR
plot_man$se_indirect = indirect_effect$se_EUR
plot_man$p_indirect = indirect_effect$p_EUR
plot_man[which(plot_man$FDR_EUR < 0.05)[2],c("Gene")] = "LRRC37A4P"


plot_man$chromosome_name <- factor(plot_man$chromosome_name, c(1:22))
colnames(plot_man)[c(20,21)] = c("CHR","BP")
# 0.05/dim(plot_man)[1] #[1] 1.006036e-05
# sig <- plot_man$p_clp_AS < 0.05/dim(plot_man)[1]

# We will use ggrepel for the annotation
library(ggrepel)

# Prepare the dataset
# Prepare the dataset
dat <- plot_man %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(plot_man, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)



nCHR <- 22

# Prepare X axis
axis.set <- dat %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#####################
## Direct


p.direct <- max(dat$p_direct[which(dat$FDR_EUR <= 0.05)]) + 0.000001
sig <- plot_man$FDR_EUR < 0.05
plot_man$Gene[which(sig == T)]
#[1] "CADM2" "LRRC37A4P"

label <- c("CADM2","LRRC37A4P")

labels_df.direct <- data.frame(label=label,
                              logP=-log10(dat$p_direct[match(label,dat$Gene)]),
                              BPcum=dat$BPcum[match(label,dat$Gene)],
                              CHR=dat$CHR[match(label,dat$Gene)])
labels_df.direct <- labels_df.direct[order(labels_df.direct$BPcum),]


manhplot.direct <- ggplot(dat, aes(x = BPcum, y = -log10(p_direct), 
                                  color = as.factor(CHR), size = -log10(p_direct))) +
  geom_point(alpha = 0.8, size=0.6) + 
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,
                     limits = c(min(dat$BPcum),max(dat$BPcum))) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 9 )) +
  scale_color_manual(values = rep(c("dark red", "red"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  geom_hline(yintercept = -log10(p.direct),
             linetype='dashed', col="black", size=0.4) +
  geom_hline(yintercept = -log10(thres),
             linetype='dashed', col="black", size=0.4) +
  guides(color = F) + 
  labs(x = NULL, 
       title = NULL) + 
  ylab( TeX("$-log_{10}(p)$") )+
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, size = 7, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 7, vjust = 0.5),
    axis.title = element_text(size=7),
    plot.title = element_text(size = 7, face = "bold"),
    plot.subtitle = element_text(size = 7)
  )+
  ggrepel::geom_label_repel(data = labels_df.direct[1,],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 1.8, segment.size = 0.2,
                            point.padding = 0.3, 
                            ylim = c(6, 10),
                            nudge_y=0.2,
                            direction = "y",segment.linetype = 1,
                            segment.curvature = -1e-20,
                            arrow = arrow(length = unit(0.015, "npc")),
                            min.segment.length = 0, force = 2,
                            box.padding = 0.5)  + 
  ggrepel::geom_label_repel(data = labels_df.direct[2,],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 1.8, segment.size = 0.2,
                            point.padding = 0.3, 
                            ylim = c(6, 10),
                            nudge_y=0.2,
                            direction = "y",segment.linetype = 1,
                            segment.curvature = -1e-20,
                            arrow = arrow(length = unit(0.015, "npc")),
                            min.segment.length = 0, force = 2,
                            box.padding = 0.5) + 
  My_Theme


#####################
## indirect

p.indirect <- p.direct

manhplot.indirect <- ggplot(dat, aes(x = BPcum, y = -log10(p_indirect), 
                                color = as.factor(CHR), size = -log10(p_indirect))) +
  geom_point(alpha = 0.8, size=0.6) + 
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,
                     limits = c(min(dat$BPcum),max(dat$BPcum))) +
  scale_y_reverse(expand = c(0,0), limits = c(9, 0))+
  scale_color_manual(values = rep(c("blue","steel blue"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  geom_hline(yintercept = -log10(p.indirect),
             linetype='dashed', col="black", size=0.4) +
  geom_hline(yintercept = -log10(thres),
             linetype='dashed', col="black", size=0.4) +
  guides(color = F) + 
  labs(x = NULL, 
       title = NULL) + 
  ylab( TeX("$-log_{10}(p)$") )+
  theme_minimal() +
 theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, size = 7, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 7, vjust = 0.5),
    axis.title = element_text(size=7),
    plot.title = element_text(size = 7, face = "bold"),
    plot.subtitle = element_text(size = 7)
  ) +
  My_Theme



p3 <- cowplot::plot_grid(manhplot.direct, manhplot.indirect, ncol=1, align="v")



###############################################################
###############################################################
###############################################################

myColors  = c("steel blue","dark red","yellow")
tmp_plot = dat[1:10,]
tmp_plot$Effect = rep(c("DE","\u03B4-IDE"),5)
tmp_plot$Effect = factor(tmp_plot$Effect,levels = c("DE","\u03B4-IDE"))
## Ancestry group color legends
tmp <- ggplot(
  data = tmp_plot,
  aes(x = Effect, y = beta_direct, ymin = se_direct, ymax = se_direct)
) +ylim(0.5,2.5)+
  geom_pointrange(aes(col = Effect), size = 0.3) +
  geom_hline(yintercept = 1, colour = "black", linetype='dashed', size=0.5) +
  xlab(NULL) +
  ylab("Relative Risk (95% Confidence Interval)") +
  geom_errorbar(aes(ymin = se_direct, ymax = se_direct, col = Effect), width = 0.2, cex = 0.2) +
  scale_color_manual(name = NULL,values = rep(c("dark red","steel blue"), 2))+
  theme_minimal() +
  theme(
    legend.key.size = unit(0.2, "mm"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text=element_text(size=7)
  )+ 
  guides(color=guide_legend(ncol = 1))


p4 <- as_ggplot(cowplot::get_plot_component(tmp, 'guide-box-right', return_all = TRUE))



p <- ggarrange(ggarrange(p1, p2,  ncol = 2, labels = c("a", "b"),font.label = list(size = 7),
                         widths = c(0.4,0.6)),
               ggarrange(p3,p4,ncol = 2, labels = c("c",NA),font.label = list(size = 7),widths = c(0.9,0.1)),
               nrow = 2,
               labels = c(NA, NA),
               heights = c(0.55,0.45)
)

ggsave(filename=paste0("Figure4.pdf"),
       plot=p, 
       path="./family/manuscript/draft4/plots/main/",
       width=180, height=160, units="mm", dpi=1000,device = cairo_pdf)

ggsave(filename=paste0("Figure4.png"),
       plot=p, device="png",
       path="./family/manuscript/draft4/plots/main/",
       width=180, height=160, units="mm", dpi=1000)

