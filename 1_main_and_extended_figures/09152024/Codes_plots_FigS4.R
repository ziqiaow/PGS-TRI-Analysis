#----------------------------------------------------------
#Figure S4
#June 2, 2024
#----------------------------------------------------------
setwd("./family/manuscript/plots/final/")
path<-"./family/manuscript/plots/final/"


library(viridis)
library(sf)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggpubr)

pgs_parents = fread("./family/simulation/UKB/plink/ea_score_parents_pruned.sscore")
load("./family/simulation/UKB/matched_parents_pheno.RData")
load("./family/simulation/UKB/selected_final.RData")
pheno = rbind(pheno_father,pheno_mother)
pheno$pgs_parents = pgs_parents$SCORE1_SUM[match(as.character(pheno$ID),as.character(pgs_parents$IID))]
pheno$bmi = selected_final$BMI[match(as.character(pheno$ID),as.character(selected_final$ID))]
pheno$townsend_index = scale(pheno$townsend_index)
pheno$education = scale(pheno$education)
pheno$bmi = scale(pheno$bmi)


#k-means clustering of geographical regions
tmp = pheno[,2:3]
tmp = scale(tmp)
set.seed(05062024)
km <- kmeans(tmp, 100)
table(km$cluster)

pheno$cluster = km$cluster

pheno.plot=pheno[,c(1,2,3,7,8,11:23)]
pheno.plot[,c(6:15)] = apply(pheno.plot[,c(6:15)],2,scale)
###################################################################################
#we make the values of each subject in the same cluster to be the same for better visualization
pheno.plot.map <- pheno.plot %>%
  summarize(PGS = mean(pgs_parents, na.rm = TRUE),
            BMI = mean(bmi, na.rm =T),
            Education = mean(education, na.rm = T),
            Townsend = mean(townsend_index, na.rm = T),
            PC1 = mean(PC1, na.rm =T),
            PC2 = mean(PC2, na.rm =T),
            PC3 = mean(PC3, na.rm =T),
            PC4 = mean(PC4, na.rm =T),
            PC5 = mean(PC5, na.rm =T),
            .by = cluster)

pheno.plot$BMI = pheno.plot.map$BMI[match(pheno.plot$cluster,pheno.plot.map$cluster)]
pheno.plot$Education = pheno.plot.map$Education[match(pheno.plot$cluster,pheno.plot.map$cluster)]
pheno.plot$PGS = pheno.plot.map$PGS[match(pheno.plot$cluster,pheno.plot.map$cluster)]
pheno.plot$Townsend = pheno.plot.map$Townsend[match(pheno.plot$cluster,pheno.plot.map$cluster)]
pheno.plot$PC1 = pheno.plot.map$PC1[match(pheno.plot$cluster,pheno.plot.map$cluster)]
pheno.plot$PC2 = pheno.plot.map$PC2[match(pheno.plot$cluster,pheno.plot.map$cluster)]
pheno.plot$PC3 = pheno.plot.map$PC3[match(pheno.plot$cluster,pheno.plot.map$cluster)]
pheno.plot$PC4 = pheno.plot.map$PC4[match(pheno.plot$cluster,pheno.plot.map$cluster)]
pheno.plot$PC5 = pheno.plot.map$PC5[match(pheno.plot$cluster,pheno.plot.map$cluster)]



# Convert to spatial points data frame
pcd_sp <- pheno.plot %>% # For object of postcodes
  st_as_sf(coords = c("birth.x", "birth.y")) %>% # Define as spatial object and identify which columns tell us the position of points
  st_set_crs(27700) # Set CRS


p1 <- ggplot(pcd_sp) + 
  geom_sf(aes(color = PGS))+
  scale_color_viridis_c(option = "A")

p2 <- ggplot(pcd_sp) + 
  geom_sf(aes(color = BMI))+
  scale_color_viridis_c(option = "A")


p3 <- ggplot(pcd_sp) + 
  geom_sf(aes(color = Education))+
  scale_color_viridis_c(option = "A")


p4 <- ggplot(pcd_sp) + 
  geom_sf(aes(color = Townsend))+
  scale_color_viridis_c(option = "A")

p_pc1 <- ggplot(pcd_sp) + 
  geom_sf(aes(color = PC1))+
  scale_color_viridis_c(option = "D")

p_pc2 <- ggplot(pcd_sp) + 
  geom_sf(aes(color = PC2))+
  scale_color_viridis_c(option = "D")

p_pc3 <- ggplot(pcd_sp) + 
  geom_sf(aes(color = PC3))+
  scale_color_viridis_c(option = "D")

p_pc4 <- ggplot(pcd_sp) + 
  geom_sf(aes(color = PC4))+
  scale_color_viridis_c(option = "D")

p_pc5 <- ggplot(pcd_sp) + 
  geom_sf(aes(color = PC5))+
  scale_color_viridis_c(option = "D")





pg1 <- ggplotGrob(p1)
pg2 <- ggplotGrob(p2)
pg3 <- ggplotGrob(p3)
pg4 <- ggplotGrob(p4)

pg_pc1 <- ggplotGrob(p_pc1)
pg_pc2 <- ggplotGrob(p_pc2)
pg_pc3 <- ggplotGrob(p_pc3)
pg_pc4 <- ggplotGrob(p_pc4)

#maxHeight = grid::unit.pmax(pg1$heights, pg2$heights)
#pg1$heights <- as.list(maxHeight)
#pg2$heights[c(7,9,11)] <- pg1$heights[c(7,9,11)] 
pg2$heights <- pg1$heights
pg3$heights <- pg1$heights
pg4$heights <- pg1$heights
pg_pc1$heights <- pg1$heights
pg_pc2$heights <- pg1$heights
pg_pc3$heights <- pg1$heights
pg_pc4$heights <- pg1$heights



p <- ggarrange(ggarrange(pg1,pg2,pg3,pg4,
                         ncol = 4, labels = c("a", "b", "c","d"),
                         widths = c(0.25,0.25,0.25,0.25)) , 
               ggarrange(pg_pc1,pg_pc2,pg_pc3,pg_pc4,
                         ncol = 4, labels = c("e", "f", "g","h"),
                         widths = c(0.25,0.25,0.25,0.25)) ,
               nrow = 2, 
               labels = c(NA, NA),
               heights = c(0.5,0.5))

ggsave(filename=paste0("Extended_Figure4_color3.pdf"),
       plot=p, 
       path="./family/manuscript/plots/final/",
       width=360, height=300, units="mm", dpi=320,device = cairo_pdf)

ggsave(filename=paste0("Extended_Figure4_color3.png"),
       plot=p, 
       path="./family/manuscript/plots/final/",
       width=360, height=300, units="mm")

