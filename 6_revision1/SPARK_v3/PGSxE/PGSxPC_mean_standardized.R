#Interaction between PGS and genetic distance based on PC (PC eigenvectors from 1kg+HGDP)
#09/03/2025
#------------------------------------------------------------------------------------------
setwd("~/work/family/data/autism/v3/GxE")

source("~/work/family/R/PGS.TRI/R/PGS-TRI.R")
source("~/work/family/R/PGS.TRI/R/pTDT.R")

#load data and data cleaning
library(data.table)
library(fastDummies)
library(ggplot2)
library(ggsci)
library(dplyr)

load("~/work/family/data/autism/v3/GxE/PRS_SPARK_PCstandardized_meanonly.RData")
load("~/work/family/data/autism/v3/roles_id_complete_casetrio_ancestry.RData")

#Define genetic distance (GD) based on top 10 PCs, use center based on EUR ASD cases in children
#For each family, the GD is fixed within family, it is an offspring exposure variable (how far away from the center of EUR ASD children)
#Original article: https://www.nature.com/articles/s41586-023-06079-4

#check reference data PC mean
load("~/work/family/data/1000G_HGDP/PCA/pc_1KG_race.RData")

#make a plot for independent ASD cases
plot_c = dat[match(id_final$subject_sp_id,dat$IID),c(1,7,11:20)]
# id = id_final$subject_sp_id[which(id_final$ancestry ==  id_final$ancestry_f  & id_final$ancestry_m == id_final$ancestry)]
# plot_c = plot_c[match(id,plot_c$IID),]
plot_c = plot_c[which(plot_c$superclass %in% c("AFR", "AMR", "EUR", "EAS", "SAS")),]

p1 <- ggplot(plot_c, aes(x = PC1, y = PC2, color = factor(superclass))) + 
  
  geom_point() +
  scale_color_jco() + theme_bw()
ggsave(filename=paste0("PCA_SPARK_ASD_children.pdf"),
       plot=p1, 
       path="~/work/family/data/1000G_HGDP/PCA/",
       width=380, height=300, units="mm", dpi=320)


#plot together with 1000G+HGDP 

#combine both dataset together
pc_1KG_race$Cohort = "1000G+HGDP"
colnames(pc_1KG_race)[16] = "superclass"
plot_c$Cohort = "SPARK"

pc_1KG_spark = rbind(pc_1KG_race[,c("IID","PC1","PC2","superclass","Cohort")],plot_c[,c("IID","PC1","PC2","superclass","Cohort")])

p1 <- ggplot(pc_1KG_spark, aes(x = PC1, y = PC2, color = factor(superclass), shape = factor(Cohort), size = factor(Cohort))) + 
  
  geom_point() +
  scale_color_jco() + 
  
  scale_shape_manual(values = c(2,20)) + 
  
  scale_size_manual(values = c(3.5,1)) + theme_bw()

ggsave(filename=paste0("PCA_SPARK_1KG_ASD_children.pdf"),
       plot=p1, 
       path="~/work/family/data/1000G_HGDP/PCA/",
       width=380, height=300, units="mm", dpi=320)


#load Finnish people
fin_ID = fread("~/work/family/data/1000G_HGDP/PCA/sample_fin.tsv")
id = intersect(fin_ID$`Sample name`,pc_1KG_race$IID)
pc_1KG_race = pc_1KG_race[match(id,pc_1KG_race$IID),]
pc_center = apply(pc_1KG_race[,c(5:14)],2,mean)

#combine both dataset together
colnames(pc_1KG_race)[16] = "superclass"
pc_1KG_race$superclass = "Finnish EUR (1000G+HGDP)"
pc_1KG_spark = rbind(pc_1KG_race[,c("IID","PC1","PC2","superclass","Cohort")],plot_c[,c("IID","PC1","PC2","superclass","Cohort")])
pc_1KG_spark$superclass = factor(pc_1KG_spark$superclass,levels = c("EUR","AFR","AMR","EAS","SAS","Finnish EUR (1000G+HGDP)"))


#Next calculate Euclidean distance of each offspring from the PC centers
GD_c = dat[match(id_final$subject_sp_id,dat$IID),c(11:20)]

#tmp1 = (GD_c)^2 #if standardized data

#GD defined by top 5 PCs
tmp1 = t(apply(GD_c, 1, function(x) (x-pc_center)^2))
GD_tmp = apply(tmp1[,c(1:5)],1,function(x) sqrt(sum(x)))
id_final$GD = GD_tmp

#GD defined by top 2 PCs
id_final$GD = GD_tmp
GD_tmp = apply(tmp1[,c(1:2)],1,function(x) sqrt(sum(x)))
id_final$GD_2pc = GD_tmp

#GD defined by top 10 PCs
GD_tmp = apply(tmp1[,c(1:10)],1,function(x) sqrt(sum(x)))
id_final$GD_10pc = GD_tmp



id_final$ancestry_test = NA
id_final$ancestry_test[which(id_final$ancestry ==  id_final$ancestry_f  & id_final$ancestry_m == id_final$ancestry)] = id_final$ancestry[which(id_final$ancestry ==  id_final$ancestry_f  & id_final$ancestry_m == id_final$ancestry)]
id_final %>%
  group_by(ancestry_test) %>%
  summarise(disp = mean(GD), sd = sd(GD))

#raw PGS
id_final$pgs_child = dat$SCORE1_SUM[match(id_final$subject_sp_id,dat$IID)]
id_final$pgs_mother = dat$SCORE1_SUM[match(id_final$biomother_sp_id ,dat$IID)]
id_final$pgs_father = dat$SCORE1_SUM[match(id_final$biofather_sp_id,dat$IID)]

#PC PGS
id_final$pgs_child_PC = dat$ASD_PGS_adj_mean_only[match(id_final$subject_sp_id,dat$IID)]
id_final$pgs_mother_PC = dat$ASD_PGS_adj_mean_only[match(id_final$biomother_sp_id ,dat$IID)]
id_final$pgs_father_PC = dat$ASD_PGS_adj_mean_only[match(id_final$biofather_sp_id,dat$IID)]

save(id_final, file="genetic_distance_meanPCstandardize_only.RData")


data.model1 = id_final[which(id_final$ancestry %in% c("EUR","AFR", "AMR",  "EAS", "SAS")),]
#18383 ASD children
table(data.model1$ancestry)

data.model1$GD = data.model1$GD / sd(data.model1$GD )

#PC mean standardized only PRS
res_int_PC=PGS.TRI(pgs_offspring = data.model1$pgs_child_PC, 
                   pgs_mother = data.model1$pgs_mother_PC, 
                   pgs_father = data.model1$pgs_father_PC,
                   GxE_int = TRUE,
                   formula = as.formula(paste("~","GD")),
                   E = data.model1, 
                   side = 2)$res_beta



data.model1$GD_2pc = data.model1$GD_2pc / sd(data.model1$GD_2pc )
PGS.TRI(pgs_offspring = data.model1$pgs_child_PC, 
        pgs_mother = data.model1$pgs_mother_PC, 
        pgs_father = data.model1$pgs_father_PC,
        GxE_int = TRUE,
        formula = as.formula(paste("~","GD_2pc")),
        E = data.model1, 
        side = 2)$Coefficients_direct


data.model1$GD_10pc = data.model1$GD_10pc / sd(data.model1$GD_10pc )
PGS.TRI(pgs_offspring = data.model1$pgs_child_PC, 
        pgs_mother = data.model1$pgs_mother_PC, 
        pgs_father = data.model1$pgs_father_PC,
        GxE_int = TRUE,
        formula = as.formula(paste("~","GD_10pc")),
        E = data.model1, 
        side = 2)$Coefficients_direct


#raw PRS
res_int_PC=PGS.TRI(pgs_offspring = data.model1$pgs_child, 
                   pgs_mother = data.model1$pgs_mother, 
                   pgs_father = data.model1$pgs_father,
                   GxE_int = TRUE,
                   formula = as.formula(paste("~","GD")),
                   E = data.model1, 
                   side = 2)$Coefficients_direct


PGS.TRI(pgs_offspring = data.model1$pgs_child, 
        pgs_mother = data.model1$pgs_mother, 
        pgs_father = data.model1$pgs_father,
        GxE_int = TRUE,
        formula = as.formula(paste("~","GD_2pc")),
        E = data.model1, 
        side = 2)$Coefficients_direct


PGS.TRI(pgs_offspring = data.model1$pgs_child, 
        pgs_mother = data.model1$pgs_mother, 
        pgs_father = data.model1$pgs_father,
        GxE_int = TRUE,
        formula = as.formula(paste("~","GD_10pc")),
        E = data.model1, 
        side = 2)$Coefficients_direct




#Make a figure for GD and RR
load("~/work/family/data/autism/v3/GxE/genetic_distance_meanPCstandardize_only.RData")
data.model1 = id_final[which(id_final$ancestry %in% c("EUR","AFR", "AMR",  "EAS", "SAS")),]
#18383 ASD children
table(data.model1$ancestry)

data.model1$GD = data.model1$GD / sd(data.model1$GD )

#group by numbers to 10 groups
breaks <- seq(min(data.model1$GD), max(data.model1$GD), length.out = 11) 
groups <- cut(data.model1$GD, breaks = breaks, labels = 1:10, include.lowest = TRUE)
data.model1$GD_group = groups
distance_by_number <- data.model1 %>%
  group_by(GD_group) %>%
  summarise(disp = mean(GD), sd = sd(GD))


ASD_effect=list()
res_direct_PC_mean = res_indirect_PC_mean= data.frame(array(0,c(10,4)))

results_all = list()
for(i in 1:10){
  data.model = data.model1[which(data.model1$GD_group == distance_by_number$GD_group[i]),]
  
  ASD_effect[[i]]=PGS.TRI(pgs_offspring = data.model$pgs_child_PC,
                          pgs_mother = data.model$pgs_mother_PC,
                          pgs_father = data.model$pgs_father_PC, 
                          GxE_int = FALSE, smalltriosize = F, parental_indirect = F)
  res_direct_PC_mean[i,] = ASD_effect[[i]]$Coefficients_direct

}
colnames(res_direct_PC_mean) = colnames(ASD_effect[[1]]$Coefficients_direct)
results_all[[1]] = res_direct_PC_mean

ASD_effect=list()
res_direct_PC_mean = res_indirect_PC_mean= data.frame(array(0,c(10,4)))

for(i in 1:10){
  data.model = data.model1[which(data.model1$GD_group == distance_by_number$GD_group[i]),]
  
  pgs_c = data.model$pgs_child_PC/sd(data.model$pgs_child_PC)
  pgs_m = data.model$pgs_mother_PC/sd(data.model$pgs_child_PC)
  pgs_f = data.model$pgs_father_PC/sd(data.model$pgs_child_PC)
  
  ASD_effect[[i]]=PGS.TRI(pgs_offspring = pgs_c,
                          pgs_mother = pgs_m,
                          pgs_father = pgs_f, 
                          GxE_int = FALSE, smalltriosize = F, parental_indirect = F)
  res_direct_PC_mean[i,] = ASD_effect[[i]]$Coefficients_direct

  
}
colnames(res_direct_PC_mean) = colnames(ASD_effect[[1]]$Coefficients_direct)
results_all[[2]] = res_direct_PC_mean



#do an interaction using this group-specific SD standardized values
data.model1$pgs_child_groupSD = 0
data.model1$pgs_mother_groupSD = 0
data.model1$pgs_father_groupSD = 0
for(i in 1:10){
  data.model = data.model1[which(data.model1$GD_group == distance_by_number$GD_group[i]),]
  data.model1$pgs_child_groupSD[which(data.model1$GD_group == distance_by_number$GD_group[i])] = data.model1$pgs_child_PC[which(data.model1$GD_group == distance_by_number$GD_group[i])]/sd(data.model$pgs_child_PC)
  data.model1$pgs_mother_groupSD[which(data.model1$GD_group == distance_by_number$GD_group[i])] = data.model1$pgs_mother_PC[which(data.model1$GD_group == distance_by_number$GD_group[i])]/sd(data.model$pgs_child_PC)
  data.model1$pgs_father_groupSD[which(data.model1$GD_group == distance_by_number$GD_group[i])] = data.model1$pgs_father_PC[which(data.model1$GD_group == distance_by_number$GD_group[i])]/sd(data.model$pgs_child_PC)
}


PGS.TRI(pgs_offspring = data.model1$pgs_child_groupSD, 
        pgs_mother = data.model1$pgs_mother_groupSD, 
        pgs_father = data.model1$pgs_father_groupSD,
        GxE_int = TRUE,
        formula = as.formula(paste("~","GD")),
        E = data.model1, 
        side = 2)$Coefficients_direct


#do an interaction using this group-specific SD standardized values and additionally assign the average value of GD to each group to correspond to Fig5
data.model1$gd_mean_by_group = 0
for(i in 1:10){
  data.model = data.model1[which(data.model1$GD_group == distance_by_number$GD_group[i]),]
  data.model1$gd_mean_by_group[which(data.model1$GD_group == distance_by_number$GD_group[i])] = mean(data.model$GD)
}


PGS.TRI(pgs_offspring = data.model1$pgs_child_groupSD, 
        pgs_mother = data.model1$pgs_mother_groupSD, 
        pgs_father = data.model1$pgs_father_groupSD,
        GxE_int = TRUE,
        formula = as.formula(paste("~","gd_mean_by_group")),
        E = data.model1, 
        side = 2)$Coefficients_direct


#add information about ancestry groups
group_summary <- data.model1 %>%
  group_by(GD_group, ancestry) %>%
  tally() %>%
  arrange(GD_group, desc(n)) %>%
  group_by(GD_group) %>%
  summarise(GD_group_ancestry = paste0(ancestry, " (N = ", n, ")", collapse = ", "))

tmp <- data.model1 %>%
  left_join(group_summary, by = "GD_group")      

table(tmp$GD_group_ancestry,tmp$GD_group)
identical(tmp$subject_sp_id,data.model1$subject_sp_id)
#[1] TRUE
data.model1$GD_group_ancestry = tmp$GD_group_ancestry


table_all = list()
for(i in 1:4){
  table_asd = results_all[[i]]
  table_asd$GD = distance_by_number$disp
  table_asd$group = distance_by_number$GD_group
  
  
  
  table_asd$OR=exp(table_asd$Estimate)
  table_asd$ci_high=exp(table_asd$Estimate+qnorm(0.975)*table_asd$Std.Error)
  table_asd$ci_low=exp(table_asd$Estimate-qnorm(0.975)*table_asd$Std.Error)
  
  
  table_plot = table_asd
  table_plot$GD_group_ancestry <- unique(data.model1$GD_group_ancestry)[order(unique(data.model1$GD_group))]
  table_plot$GD_group_ancestry <- factor(table_plot$GD_group_ancestry,levels = unique(table_plot$GD_group_ancestry))

  table_plot$sig = "No"
  table_plot$sig[which(table_plot$Pvalue < 0.05)] = "Yes"
  table_plot$sig = factor(table_plot$sig,levels = c("Yes","No"))
  
  table_all[[i]] = table_plot
}

names(table_all)=c("PCmeanstandardizeonly","PCmeanstandardizeonly_SDbyGDgroup","raw","raw_SDbyGDgroup")
save(table_all,file="PGS_GD_plot_PCmeanstandardizeonly_raw.RData")
table_output = do.call(rbind,table_all)
write.csv(table_output,file="PGS_GD_plot_PCmeanstandardizeonly_raw.csv")





#Plot PGSxPC interactions
load("~/work/family/data/autism/v3/GxE/PGS_GD_plot_PCmeanstandardizeonly_raw.RData")

# Reduce the opacity of the grid lines: Default is 255
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)

p3 <- ggplot(
  data = table_all[[2]],
  aes(x = GD, y = OR, ymin = ci_low, ymax = ci_high))+  
  geom_pointrange(aes(col = GD_group_ancestry,shape = sig), position=position_dodge(width=0.8)) +scale_shape_manual(values = c(8, 16),name="P < 0.05")+
  geom_hline(yintercept = 1, colour = "black", linetype='dashed', size=0.8) +
  xlab("Mean Genetic Distance from Reference Dataset") +
  ylab("Relative Risk for ASD (95% CI)") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, col = GD_group_ancestry), width = 0.2, cex = 0.8, position=position_dodge(width=0.5)) +
  scale_color_ucscgb() + theme_classic() +
  theme(
    panel.background = element_blank(), strip.background = element_rect(colour = NA, fill = NA),
    strip.text.y = element_text(face = "bold",size = 10, color = "black"),
    panel.grid.major.y = element_line(colour = col_grid, size = 1 ),
    strip.text = element_text(face = "bold",size = 10, color = "black"),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "right",
    axis.text = element_text(face = "bold",size = 10, color = "black"),
    axis.title = element_text(face = "bold",size = 10, color = "black",hjust = 0.5),
    plot.title = element_text(face = "bold", hjust = 0.5,size = 10, color = "black"),
    legend.title=element_text(face = "bold", size=10), 
    legend.text=element_text(face = "bold", size=10)
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
    size = 3        
  )


ggsave(filename=paste0("RR_GD_and_ancestry_PCmeanstandardizeonly_SDbyGDgroup.png"),
       plot=p2, device="png",
       path="~/work/family/manuscript/plots/draft2/",
       width=300, height=150, units="mm", dpi=320)


