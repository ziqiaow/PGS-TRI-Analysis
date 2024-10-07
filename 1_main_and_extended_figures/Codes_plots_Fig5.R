
## Fig 5

library(patchwork)
library(dplyr, warn = FALSE)
library(ggplot2)
library(RColorBrewer)

# Panel A

load("./family/data/oral_cleft/results/interaction_08292024.RData")
load("./family/data/oral_cleft/results/nurture_08292024.RData")
load("./family/data/oral_cleft/results/PGS_main_cpt_ptdt_08292024.RData")


#plot only CL/P, CL, and CL&P subtype
colnames(res_interaction)[7:9] = colnames(pgscpt_main)[7:9]
pgscpt_main = data.frame(pgscpt_main)
pgscpt_nurture = data.frame(pgscpt_nurture)
pgscpt_main$subtype = pgscpt_nurture$subtype = c("All","CL/P","CL","CL&P","CP")
pgscpt_main$envir = "DE"
pgscpt_nurture$envir = "\u03B4-IDE"
res_interaction$subtype[which(res_interaction$subtype == "all")]="All"
res_interaction$envir[which(res_interaction$envir == "smoke")]="PGS x Smoking"
res_interaction$envir[which(res_interaction$envir == "ets")]="PGS x Envir Smoke"
res_interaction$envir[which(res_interaction$envir == "sex")]="PGS x Sex (Female)"

table_clp=rbind(pgscpt_main[c(2,3,4),],res_interaction[-which(res_interaction$subtype=="CP" | res_interaction$subtype=="All" |res_interaction$envir =="vit" | res_interaction$envir =="drink"),],pgscpt_nurture[c(2,3,4),]) #Remove smoke, vitamin from the figure

table_clp=unname(table_clp)
table_clp2=rbind(data.frame(table_clp[,1:3]),data.frame(table_clp[,4:6]),data.frame(table_clp[,7:9]))
table_clp2$Race=c(rep("Cross Population",15),rep("EUR",15),rep("Asian",15))
table_clp2$Subtype = rep(table_clp[,10],3)
table_clp2$Effect = rep(table_clp[,11],3)

table_clp2$OR=exp(table_clp2$X1)
table_clp2$ci_high=exp(table_clp2$X1+qnorm(0.975)*table_clp2$X2)
table_clp2$ci_low=exp(table_clp2$X1-qnorm(0.975)*table_clp2$X2)

#main effect
table_plot = table_clp2
effect_order = c("DE","\u03B4-IDE","PGS x Smoking", "PGS x Envir Smoke","PGS x Sex (Female)")
table_plot$Effect = factor(table_plot$Effect, levels = rev(effect_order))
table_plot$Subtype = factor(table_plot$Subtype, levels = rev(c("CL/P","CL","CL&P")))
table_plot$Race = factor(table_plot$Race, levels = c("Cross Population","EUR","Asian"))
table_plot$sig = "No"
table_plot$sig[which(table_plot$X3 < 0.05)] = "Yes"

# Reduce the opacity of the grid lines: Default is 255
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

p1<-ggplot(
  data = table_plot,
  aes(x = Effect, y = OR, ymin = ci_low, ymax = ci_high, group = Subtype)
) +
  geom_pointrange(aes(col = Subtype,shape = sig), position=position_dodge(width=0.6)) +scale_shape_manual(values = c(16, 8),name="P < 0.05")+
  geom_hline(yintercept = 1, colour = "grey", linetype='dashed', size=0.5) +
  xlab(NULL) +
  ylab("Relative Risk (95% CI)") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, col = Subtype), width = 0.25, cex = 0.7, position=position_dodge(width=0.6)) +
  scale_color_manual(values = c("green",brewer.pal(9,"YlGnBu")[c(5,7)],"dark green"))+#rep(c("steel blue","dark red","orange","purple"), 1)) + #"orange","purple"
  theme_classic() +
  facet_grid(~Race)+#,  scales = "free") +
  theme(
    panel.background = element_blank(), strip.background = element_rect(colour = NA, fill = NA),
    strip.text.y = element_text(face = "bold",size = 10, color = "black"),
    panel.grid.major.y = element_line(colour = col_grid, size = 1 ),
    strip.text = element_text(face = "bold",size = 10, color = "black"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "right",
    axis.text = element_text(face = "bold",size = 10, color = "black"),
    axis.title = element_text(face = "bold",size = 10, color = "black",hjust = 0.5),
    plot.title = element_text(face = "bold", hjust = 0.5,size = 10, color = "black")
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
load("./phewas/plink/oc/summary/plots/RNAseq_maineffect_anno_redo.RData")

plot_man=prot.anno
plot_man$chromosome_name <- factor(plot_man$chromosome_name, c(1:22))
colnames(plot_man)[c(28,29)] = c("CHR","BP")

# We will use ggrepel for the annotation
library(ggrepel)

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
## Asian


p.asian <- 0.05/dim(plot_man)[1]
sig <- plot_man$p_clp_AS < 0.05/dim(plot_man)[1]
plot_man$Gene[which(sig == T)]

label <- c("TRAF3IP3")

labels_df.asian <- data.frame(label=label,
                             logP=-log10(dat$p_clp_AS[match(label,dat$Gene)]),
                             BPcum=dat$BPcum[match(label,dat$Gene)],
                             CHR=dat$CHR[match(label,dat$Gene)])
labels_df.asian <- labels_df.asian[order(labels_df.asian$BPcum),]


manhplot.asian <- ggplot(dat, aes(x = BPcum, y = -log10(p_clp_AS), 
                                 color = as.factor(CHR), size = -log10(p_clp_AS))) +
  geom_point(alpha = 0.8, size=0.8) + 
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,
                     limits = c(min(dat$BPcum),max(dat$BPcum))) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 12 )) +
  scale_color_manual(values = rep(c("#4292c6", "#08306b"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  geom_hline(yintercept = -log10(p.asian),
             linetype='dashed', col="red", size=0.5) +
  geom_hline(yintercept = -log10(thres),
             linetype='dashed', col="blue", size=0.5) +
  guides(color = F) + 
  labs(x = NULL, 
       title = NULL) + 
  ylab( TeX("$-log_{10}(p)$") )+
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, size = 5, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 6, vjust = 0.5),
    axis.title = element_text(size=7),
    plot.title = element_text(size = 7, face = "bold"),
    plot.subtitle = element_text(size = 7)
  ) + 
  ggrepel::geom_label_repel(data = labels_df.asian[1,],
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col="black",
                            size = 3, segment.size = 0.5,
                            point.padding = 0.5, 
                            ylim = c(6, 10),
                            nudge_y=0.2,
                            direction = "y",segment.linetype = 1,
                            segment.curvature = -1e-20,
                            arrow = arrow(length = unit(0.015, "npc")),
                            min.segment.length = 0, force = 2,
                            box.padding = 0.5) + 
  My_Theme


#####################
## EUR

p.eur <- 0.05/dim(plot_man)[1]

label <- c("TRAF3IP3")

labels_df.eur <- data.frame(label=label,
                              logP=-log10(dat$p_clp_EUR[match(label,dat$Gene)]),
                              BPcum=dat$BPcum[match(label,dat$Gene)],
                              CHR=dat$CHR[match(label,dat$Gene)])
labels_df.eur <- labels_df.eur[order(labels_df.eur$BPcum),]

manhplot.eur <- ggplot(dat, aes(x = BPcum, y = -log10(p_clp_EUR), 
                                  color = as.factor(CHR), size = -log10(p_clp_EUR))) +
  geom_point(alpha = 0.8, size=0.8) + 
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center,
                     limits = c(min(dat$BPcum),max(dat$BPcum))) +
  scale_y_reverse(expand = c(0,0), limits = c(12, 0))+
  scale_color_manual(values = rep(c("green","dark green"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  geom_hline(yintercept = -log10(p.eur),
             linetype='dashed', col="red", size=0.5) +
  geom_hline(yintercept = -log10(thres),
             linetype='dashed', col="blue", size=0.5) +
  guides(color = F) + 
  labs(x = NULL, 
       title = NULL) + 
  ylab( TeX("$-log_{10}(p)$") )+
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 6, vjust = 0.5),
    axis.title = element_text(size=7),
    plot.title = element_text(size = 7, face = "bold"),
    plot.subtitle = element_text(size = 7)
  )+
  ggrepel::geom_label_repel(data = labels_df.eur,
                            aes(x = .data$BPcum,
                                y = .data$logP,
                                label = .data$label), col=c("black"),
                            size = 3, segment.size = 0.5,
                            point.padding = 0.5, 
                            direction = "y",
                            ylim = c( -10, -6),
                            min.segment.length = 0, force = 5,segment.linetype = 1,
                            segment.curvature = -1e-20,
                            arrow = arrow(length = unit(0.015, "npc")),
                            nudge_y = 1.2,
                            box.padding = 0.5)+
My_Theme



p3 <- cowplot::plot_grid(manhplot.asian, manhplot.eur, ncol=1, align="v")



###############################################################
###############################################################
###############################################################

myColors  = c("dark green","#08306b","yellow")
tmp_plot = table_plot
tmp_plot = table_plot[-which(table_plot$Race == "Cross Population"),]
tmp_plot$Race = factor(tmp_plot$Race,levels = c("Asian","EUR"))
## Ancestry group color legends
tmp <- ggplot(
  data = tmp_plot,
  aes(x = Race, y = X1, ymin = ci_low, ymax = ci_high)
) +ylim(0.5,2.5)+
  geom_pointrange(aes(col = Race)) +
  geom_hline(yintercept = 1, colour = "black", linetype='dashed', size=0.5) +
  xlab(NULL) +
  ylab("Relative Risk (95% Confidence Interval)") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, col = Race), width = 0.2, cex = 0.5) +
  scale_color_manual(name = NULL,values = rep(c("#08306b","dark green"), 2))+
  theme_minimal() +
  theme(
    legend.key.size = unit(2, "mm"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text=element_text(size=10)
  )+ 
  guides(color=guide_legend(ncol = 1))


p4 <- as_ggplot(cowplot::get_plot_component(tmp, 'guide-box-right', return_all = TRUE))



###############################################################
###############################################################
###############################################################


p <- ggarrange(p1,
               ggarrange(p3,p4,ncol = 2, labels = c("b",NA),widths = c(0.9,0.1)),
               nrow = 2, 
               labels = c("a", NA),
               heights = c(0.65,0.35)
)

ggsave(filename=paste0("Figure5_09152024.pdf"),
       plot=p, 
       path="./family/manuscript/plots/final/",
       width=180, height=180, units="mm", dpi=600,device = cairo_pdf)

ggsave(filename=paste0("Figure5_09152024.png"),
       plot=p, 
       path="./family/manuscript/plots/final/",
       width=180, height=180, units="mm", dpi=600)

