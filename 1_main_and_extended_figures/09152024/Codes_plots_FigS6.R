#Extended Figure 6
#PCA analysis of the 27 nightingale markers in 1000G
setwd("./phewas/1000G")
library("ggplot2")
library("ggfortify")
library(ggpubr)
library(data.table)
type="Nightingale"
load("./phewas/1000G/Nightingale_data.RData")
scores_nightingale = data.frame(scores_nightingale)
info=fread(paste0("./phewas/data/",type,"_trait_validation_results_with_OMICSPRED_ID.csv"))
rm_name=unique(info$`OMICSPRED ID`[which(info$Internal_R2<0.1 | info$`#SNP` < 5)])
id = intersect(rm_name,colnames(scores_nightingale))
#remove those with SNPs<5
scores_nightingale=scores_nightingale[,-match(rm_name,colnames(scores_nightingale))]
rownames(scores_nightingale) = scores_nightingale$IID
scores_nightingale = scores_nightingale[,-1]
info_clean = info[match(colnames(scores_nightingale),info$`OMICSPRED ID`),1:6]

pc <- prcomp(t(scores_nightingale),
             center = TRUE,
             scale. = TRUE)
attributes(pc)
p1 <- autoplot(pc,data = info_clean,colour = "Group")+geom_point(pch=21, colour="Black")
p2 <- autoplot(pc,data = info_clean,colour = "Group", x = 1, y = 3)+geom_point(pch=21, colour="Black")
p3 <- autoplot(pc,data = info_clean,colour = "Group", x = 2, y = 3)+geom_point(pch=21, colour="Black")
p4 <- autoplot(pc,data = info_clean,colour = "Group", x = 1, y = 4)+geom_point(pch=21, colour="Black")


p1 <- p1+ theme(legend.position="none")
p2 <- p2+ theme(legend.position="none")



# compute total variance
variance = (pc$sdev^2 / sum(pc$sdev^2))
p5<- qplot(c(1:4), variance[1:4]) + 
  geom_col(col = "black",fill="dark green")+
  xlab("Principal Component") + 
  ylab("Prop of Variance Explained") +
  ggtitle(NULL) +
  ylim(0, 1) 

p <- ggarrange(ggarrange(p1, p2,
                         ncol = 2, labels = c("a", "b"),
                         widths = c(0.45,0.45)) , 
               ggarrange(p4,p5,ncol = 2, labels = c("c","d"),widths = c(0.45,0.45),common.legend = T,legend="bottom"),
               nrow = 2, 
               labels = c(NA, NA),
               heights = c(0.45,0.47)
)


ggsave(filename=paste0("Extended_Figure6.pdf"),
       plot=p, device="pdf",
       path="./family/manuscript/plots/final/",
       width=160, height=160, units="mm", dpi=320)
ggsave(filename=paste0("Extended_Figure6.png"),
       plot=p, device="png",
       path="./family/manuscript/plots/final/",
       width=160, height=160, units="mm", dpi=320)


