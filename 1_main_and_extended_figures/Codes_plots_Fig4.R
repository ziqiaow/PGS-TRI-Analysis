
## Fig 4

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

setwd("./family/data/autism/results")
load("./family/data/autism/results/PGS_main_pgscpt.RData")
load("./family/data/autism/results/nurture.RData")

main=data.frame(tmp_res_main[c(2:5,1),])
main$race=c("European","African","America",
            "EAS/SAS","Cross Population")

nurture=data.frame(pgscpt_nurture[c(2:5,1),])
nurture$race=c("European","African","America","EAS/SAS","Cross Population")

table_asd=main

table_asd$effect = c("DE")



table_asd$OR=exp(table_asd$beta)
table_asd$ci_high=exp(table_asd$beta+qnorm(0.975)*table_asd$se)
table_asd$ci_low=exp(table_asd$beta-qnorm(0.975)*table_asd$se)


table_plot = table_asd
table_plot$race = factor(table_plot$race, levels = rev(unique(table_plot$race)))

table_plot1= table_plot

table_asd=nurture

table_asd$effect = c("\u03B4-IDE")




table_asd$OR=exp(table_asd$beta)
table_asd$ci_high=exp(table_asd$beta+qnorm(0.975)*table_asd$se)
table_asd$ci_low=exp(table_asd$beta-qnorm(0.975)*table_asd$se)


table_plot = table_asd
table_plot$race = factor(table_plot$race, levels = rev(unique(table_plot$race)))

table_plot_final = rbind(table_plot,table_plot1)
table_plot_final$sig = "no"
table_plot_final$sig[which(table_plot_final$p < 0.05)] = "yes"
table_plot_final$effect = factor(table_plot_final$effect,levels = rev(c("DE","\u03B4-IDE")))


#main effect
# Reduce the opacity of the grid lines: Default is 255
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)
table_plot_final$trait = "PGS for ASD"
p1<-ggplot(
  data = table_plot_final,
  aes(x = race, y = OR, ymin = ci_low, ymax = ci_high, group = effect)
) +ylim(0.5,2.5)+
  geom_pointrange(aes(col = effect,shape = sig), position=position_dodge(width=0.5)) +scale_shape_manual(values = c(16, 8))+
  geom_hline(yintercept = 1, colour = "black", linetype='dashed', size=0.5) +
  xlab(NULL) +
  ylab("Relative Risk for ASD (95% CI)") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, col = effect), width = 0.2, cex = 0.5, position=position_dodge(width=0.5)) +
  scale_color_manual(values = rep(c("steel blue","dark red"), 1)) + #"orange","purple"
theme_classic() +
  facet_wrap(~trait, strip.position = "top", nrow = 2, scales = "free_y") +
  theme(
    panel.background = element_blank(), strip.background = element_rect(colour = NA, fill = NA),
    strip.text.y = element_text(face = "bold",size = 8, color = "black"),
    panel.grid.major.y = element_line(colour = col_grid, size = 1 ),
    strip.text = element_text(face = "bold",size = 8, color = "black"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "none",
    axis.text = element_text(face = "bold",size = 8, color = "black"),
    axis.title = element_text(face = "bold",size = 8, color = "black",hjust = 0.5),
    plot.title = element_text(face = "bold", hjust = 0.5,size = 8, color = "black")
  ) + coord_flip()





# Panel B

setwd("./family/data/autism/results")
load("./family/data/autism/results/PGS_main_pgscpt.RData")
load("./family/data/pleiotropy/results_all_08212024.RData")

res_plot=delta_all[-c(6,8,9,12),c(16:18,4:6)]
res_plot=res_plot[c(1,2,3,5,4,7,8,6),]

colnames(res_plot)=rep(colnames(main)[1:3],2)
table_asd=cbind(res_plot[,1:3],race = "PGS for Other Traits") 

table_asd$Trait = c("Educational Attainment","Schizophrenia","Major Depression","ADHD",
                    "Bipolar Disorder",
                    "Neuroticism","Insomnia","BMI")



table_asd$OR=exp(table_asd$beta)
table_asd$ci_high=exp(table_asd$beta+qnorm(0.975)*table_asd$se)
table_asd$ci_low=exp(table_asd$beta-qnorm(0.975)*table_asd$se)

table_plot = table_asd
table_plot$Trait = factor(table_plot$Trait, levels = rev(c("Educational Attainment","Schizophrenia","Major Depression","ADHD",
                                                           "Bipolar Disorder","Neuroticism","Insomnia",
                                                           "BMI")))
table_plot$sig = "No"
table_plot$sig[which(table_plot$p < 0.05)] = "Yes"

table_plot1 = table_plot



res_plot=results_all[-c(6,8,9,12),c(16:18,4:6)]
res_plot=res_plot[c(1,2,3,5,4,7,8,6),]

colnames(res_plot)=rep(colnames(main)[1:3],2)
table_asd=cbind(res_plot[,1:3],race = "PGS for Other Traits") 
table_asd$Trait = c("Educational Attainment","Schizophrenia","Major Depression","ADHD",
                    "Bipolar Disorder",
                    "Neuroticism","Insomnia","BMI")



table_asd$OR=exp(table_asd$beta)
table_asd$ci_high=exp(table_asd$beta+qnorm(0.975)*table_asd$se)
table_asd$ci_low=exp(table_asd$beta-qnorm(0.975)*table_asd$se)

table_plot = table_asd
table_plot$Trait = factor(table_plot$Trait, levels = rev(c("Educational Attainment","Schizophrenia","Major Depression","ADHD",
                                                           "Bipolar Disorder","Neuroticism","Insomnia",
                                                           "BMI")))
table_plot$sig = "No"
table_plot$sig[which(table_plot$p < 0.05)] = "Yes"

table_plot_total = rbind(table_plot,table_plot1)
table_plot_total$effect = c(rep("DE",8),rep("\u03B4-IDE",8))
table_plot_total$effect = factor(table_plot_total$effect,levels = rev(c("DE","\u03B4-IDE")))
#main effect
# Reduce the opacity of the grid lines: Default is 255
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

p2 <- ggplot(
  data = table_plot_total,
  aes(x = Trait, y = OR, ymin = ci_low, ymax = ci_high, group = effect)
)+ylim(0.8,1.35) +
  geom_pointrange(aes(col = effect,shape = sig), position=position_dodge(width=0.5)) +scale_shape_manual(values = c(16, 8),name="P < 0.05")+
  geom_hline(yintercept = 1, colour = "black", linetype='dashed', size=0.5) +
  xlab(NULL) +
  ylab("Relative Risk for ASD (95% CI)") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, col = effect), width = 0.2, cex = 0.5, position=position_dodge(width=0.5)) +
  scale_color_manual(values = rep(c("steel blue","dark red"), 1),name="Genetic Effect") + #"orange","purple"
  facet_wrap(~race, strip.position = "top", nrow = 2, scales = "free_y") +
  theme_classic() +
  theme(
    panel.background = element_blank(), strip.background = element_rect(colour = NA, fill = NA),
    strip.text.y = element_text(face = "bold",size = 8, color = "black"),
    panel.grid.major.y = element_line(colour = col_grid, size = 1 ),
    strip.text = element_text(face = "bold",size = 8, color = "black"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "right",
    axis.text = element_text(face = "bold",size = 8, color = "black"),
    axis.title = element_text(face = "bold",size = 8, color = "black",hjust = 0.5),
    plot.title = element_text(face = "bold", hjust = 0.5,size = 8, color = "black")
  ) + coord_flip() + 
  guides(color = guide_legend(reverse = TRUE),shape = guide_legend(reverse = TRUE))


########################################################################
########################################################################
########################################################################
#Panel C
#Make a heatmap for the metabolites DNE
type="Nightingale"
load(paste0("./phewas/summary/redo_06122024/",type,"_results.RData"))
res_nurture=data.frame(res_nurture)

info=fread(paste0("./phewas/data/",type,"_trait_validation_results_with_OMICSPRED_ID.csv"))
rm_name=info$`OMICSPRED ID`[which(info$Internal_R2<0.1 | info$`#SNP`<5)]
res_nurture=res_nurture[-match(intersect(rm_name,rownames(res_nurture)),rownames(res_nurture)),]
nightingale_results_nurture=res_nurture
nightingale_results_nurture$ID=rownames(nightingale_results_nurture)
nightingale_results_nurture=merge(nightingale_results_nurture,info[,1:6],by.x = "ID",by.y="OMICSPRED ID")


nightingale_results_nurture=nightingale_results_nurture[,c(2,4,8,10,20,21,22,23)]
nightingale_results_nurture$Subgroup=gsub("\\(.*","",nightingale_results_nurture$Subgroup)
nightingale_results_nurture$Subgroup = trimws(nightingale_results_nurture$Subgroup)
tmp1=nightingale_results_nurture[,c(1,3,6,7,8)]
colnames(tmp1)=c("All","European","Metabolites","Group","Subgroup")
long <- melt(setDT(tmp1), id.vars = c("Metabolites","Group","Subgroup"), variable.name = "beta")

tmp2=nightingale_results_nurture[,c(2,4,6,7,8)]
colnames(tmp2)=c("All","European","Metabolites","Group","Subgroup")
long2 <- melt(setDT(tmp2), id.vars = c("Metabolites","Group","Subgroup"), variable.name = "p")
colnames(long)[c(4,5)]=c("population","beta")
colnames(long2)[c(4,5)]=c("population","p")
plot_metabolite=merge(long,long2,by=c("Metabolites","Group","Subgroup","population"))
plot_metabolite$population = factor(plot_metabolite$population,levels = rev(unique(plot_metabolite$population)))

plot_metabolite$Subgroup[which(plot_metabolite$Subgroup == "Very small VLDL")]="VLDL"
plot_metabolite$Subgroup[which(plot_metabolite$Subgroup == "small LDL")]="sdLDL"
plot_metabolite = plot_metabolite[order(plot_metabolite$Group,plot_metabolite$Subgroup),]
tmp_row  <- plot_metabolite[c(1,2),]
plot_metabolite[c(1,2),]  <- plot_metabolite[c(5,6),]
plot_metabolite[c(5,6),]  <- tmp_row
plot_metabolite$Metabolites <- factor(plot_metabolite$Metabolites , levels = unique(plot_metabolite$Metabolites ))


label = plot_metabolite$Group
p3 <- plot_metabolite %>% 
  ggplot(aes(x = Metabolites, y = population, fill = beta, 
             label = ifelse(p < 0.05/2, "*", "")) ) + 
  geom_tile() +  geom_text(vjust = 0.5,hjust=0.5,size=5)+
  ylab(NULL) +xlab(NULL)+
  scale_fill_viridis(option = "C",limits=c(-0.2,0.2),name="\u03B4-IDE")+
  theme_classic() +
  theme(
    panel.background = element_blank(), strip.background = element_rect(colour = NA, fill = NA),
    legend.position = "right",
    axis.text.x = element_text(face = "bold",size = 8, color = "black",angle = 45,hjust=1, vjust=1),
    axis.text.y = element_text(face = "bold",size = 8, color = "black"),
    axis.title = element_text(face = "bold",size = 8, color = "black",hjust = 0.5)
  ) + theme(line = element_blank())



default.colors <- c(pal_nejm("default")(5)[c(1:4)],cols <- brewer.pal(9,"YlGnBu")[c(4:7,9)])
names(x = default.colors) <- c(unique(plot_metabolite$Group),unique(plot_metabolite$Subgroup)[-1])

tmp2 = cbind( default.colors[plot_metabolite$Group],default.colors[plot_metabolite$Subgroup])

# color only
data <- unique(plot_metabolite$Group)
colorScales <- default.colors[1:4]
names(colorScales) <- data

h1 <- ggplot() +
  geom_point(aes(x = data, y = data, color = as.character(data)),
             size=5,shape = 15
  ) +
  scale_color_manual(
    name = "Group",
    values = colorScales
  ) +
  theme(legend.position = "bottom")


data2 <- unique(plot_metabolite$Subgroup)[-1]
colorScales <- default.colors[5:10]
names(colorScales) <- data2


h2 <- ggplot() +
  geom_point(aes(x = data2, y = data2, color = as.character(data2)),
             size=5,shape = 15
  ) +
  scale_color_manual(
    name = "Lipoprotein Subgroup",
    values = colorScales
  ) +
  theme(legend.position = "bottom")



leg1 <- cowplot::get_plot_component(h1, 'guide-box-bottom', return_all = TRUE)
leg2 <- cowplot::get_plot_component(h2, 'guide-box-bottom', return_all = TRUE)

# combine all legends
leg12 <- ggarrange(leg1,leg2, 
                    nrow = 2,labels = c("", "",""),
                    heights = c(0.5,0.5))



pbuild <- ggplot_build(plot = p3)
# scale the height of the bar
y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.005
y.max <- y.pos + 0.1 * y.range
x.min <- min(pbuild$layout$panel_params[[1]]$x.range) + 0.1
x.max <- max(pbuild$layout$panel_params[[1]]$x.range) - 0.1

p3 <- p3+annotation_raster(raster = t(x = tmp2),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max)+  coord_cartesian(ylim = c(1, 2.2), clip = "off")+
  annotation_custom(grob = grid::textGrob(label = "Group", hjust = 1.05, gp = grid::gpar(cex = 0.65)), ymin = mean(c(y.pos, y.max))+0.07, ymax = mean(c(y.pos, y.max))+0.07, xmin = -Inf, xmax = -Inf)+
  annotation_custom(grob = grid::textGrob(label = "Subgroup", hjust = 0.99, gp = grid::gpar(cex = 0.65)), ymin = mean(c(y.pos, y.max))-0.07, ymax = mean(c(y.pos, y.max))-0.07, xmin = -Inf, xmax = -Inf) 

p3 <- p3+ggtitle("ASD Risk associated with \u03B4-IDE of obesity-related metabolite PGSs")+ theme(plot.title = element_text(size = 10,face = "bold"))

final_p <- ggarrange( leg12, p3, nrow = 2,labels = c("",""),heights = c(0.05,0.9))



p <- ggarrange(ggarrange(ggarrange(p1, p2,
                                   ncol = 2, labels = c("a", "b"),
                                   widths = c(0.4,0.6)),
                         p3,leg12,nrow = 3,labels = c("","c",""),heights = c(0.5,0.45,0.1)
))

ggsave(filename=paste0("Figure4_09252024.pdf"),
       plot=p, 
       path="./family/manuscript/plots/final/",
       width=180, height=200, units="mm", dpi=320,device = cairo_pdf)

ggsave(filename=paste0("Figure4_09252024.png"),
       plot=p, device="png",
       path="./family/manuscript/plots/final/",
       width=180, height=200, units="mm", dpi=320)
