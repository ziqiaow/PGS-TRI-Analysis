#Interaction between PGS and genetic distance based on PC (PC eigenvectors from 1kg+HGDP)
#03/11/2025
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

load("~/work/family/data/autism/v3/PRS_SPARK_PCstandardized.RData")
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

tmp1 = t(apply(GD_c, 1, function(x) (x-pc_center)^2))
GD_tmp = apply(tmp1[,c(1:5)],1,function(x) sqrt(sum(x)))
id_final$GD = GD_tmp
id_final$ancestry_test = NA
id_final$ancestry_test[which(id_final$ancestry ==  id_final$ancestry_f  & id_final$ancestry_m == id_final$ancestry)] = id_final$ancestry[which(id_final$ancestry ==  id_final$ancestry_f  & id_final$ancestry_m == id_final$ancestry)]
id_final %>%
  group_by(ancestry_test) %>%
  summarise(disp = mean(GD), sd = sd(GD))

#PC PGS
id_final$pgs_child_PC = dat$ASD_PGS_adj_mean[match(id_final$subject_sp_id,dat$IID)]
id_final$pgs_mother_PC = dat$ASD_PGS_adj_mean[match(id_final$biomother_sp_id ,dat$IID)]
id_final$pgs_father_PC = dat$ASD_PGS_adj_mean[match(id_final$biofather_sp_id,dat$IID)]

#PC sd PGS
id_final$pgs_child_PC_sd = dat$ASD_PGS_adj_sd[match(id_final$subject_sp_id,dat$IID)]
id_final$pgs_mother_PC_sd = dat$ASD_PGS_adj_sd[match(id_final$biomother_sp_id ,dat$IID)]
id_final$pgs_father_PC_sd = dat$ASD_PGS_adj_sd[match(id_final$biofather_sp_id,dat$IID)]

save(id_final, file="genetic_distance.RData")


data.model1 = id_final[which(id_final$ancestry %in% c("EUR","AFR", "AMR",  "EAS", "SAS")),]
#18383 ASD children
table(data.model1$ancestry)
#  AFR   AMR   EAS   EUR   SAS 
#1235  2410   442 13668   628 
data.model1$GD = data.model1$GD / sd(data.model1$GD )
res_int_PC=PGS.TRI(pgs_offspring = data.model1$pgs_child_PC, 
                   pgs_mother = data.model1$pgs_mother_PC, 
                   pgs_father = data.model1$pgs_father_PC,
                   GxE_int = TRUE,
                   formula = as.formula(paste("~","GD")),
                   E = data.model1, 
                   side = 2)$res_beta


#export to latex table
PCs = c(2,5,10)
res_int_PC = list()
for(i in 1:3){
  tmp1 = t(apply(GD_c, 1, function(x) (x-pc_center)^2))
  GD_tmp = apply(tmp1[,c(1:PCs[i])],1,function(x) sqrt(sum(x)))
  id_final$GD = GD_tmp
  
  data.model1 = id_final[which(id_final$ancestry %in% c("EUR","AFR", "AMR",  "EAS", "SAS")),]
  #18383 ASD children
  table(data.model1$ancestry)
  #  AFR   AMR   EAS   EUR   SAS 
  #1235  2410   442 13668   628 
  data.model1$GD = data.model1$GD / sd(data.model1$GD )
  res_int_PC[[i]]=PGS.TRI(pgs_offspring = data.model1$pgs_child_PC, 
                     pgs_mother = data.model1$pgs_mother_PC, 
                     pgs_father = data.model1$pgs_father_PC,
                     GxE_int = TRUE,
                     formula = as.formula(paste("~","GD")),
                     E = data.model1, 
                     side = 2)$res_beta
  
  
  
}

output_tab = array(0,c(3,3))
for(i in 1:3){
  output_tab[i,] = res_int_PC[[i]][2,c(1,2,4)]
}
rownames(output_tab) = c(paste0("Top ",c(2,5,10)," PCs"))
colnames(output_tab) = colnames(res_int_PC[[1]])[c(1,2,4)]
library(xtable)
xtable(output_tab,digits = c(3,3,3,6))
