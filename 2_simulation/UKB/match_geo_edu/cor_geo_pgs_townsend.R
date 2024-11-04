#check the correlation between parental PGS and townsend index, edu, height, BMI by clusters of geographical regions
#load parental PGS scores to check the correlations between EA pgs and birth locations
pgs_parents = fread("plink/ea_score_parents_pruned.sscore")
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

#calculate the means for PGS and towsend index and edu by clusters (between-cluster correlation)
means = aggregate(pheno[,c("education","townsend_index","pgs_parents","bmi")], list(km$cluster), mean, na.rm = TRUE)
cor.test(means$townsend_index,means$pgs_parents)
# 	Pearson's product-moment correlation
# 
# data:  means$townsend_index and means$pgs_parents
# t = -0.90642, df = 98, p-value = 0.3669
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2825391  0.1071561
# sample estimates:
#   cor 
# -0.09118113 


cor.test(means$education,means$pgs_parents)
#	Pearson's product-moment correlation
# 
# data:  means$education and means$pgs_parents
# t = -7.0092, df = 98, p-value = 3.076e-10
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.6953507 -0.4302744
# sample estimates:
#   cor 
# -0.5778558 

cor.test(means$bmi,means$pgs_parents)
#		Pearson's product-moment correlation
# 
# data:  means$bmi and means$pgs_parents
# t = -5.3765, df = 98, p-value = 5.153e-07
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.6159423 -0.3098975
# sample estimates:
#   cor 
# -0.4772647 

cor.test(means$bmi,means$education)
#	Pearson's product-moment correlation
# 
# data:  means$bmi and means$education
# t = 3.7667, df = 98, p-value = 0.0002824
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.1711586 0.5159967
# sample estimates:
#   cor 
# 0.3556212 






#within cluster correlation
pheno$km_cluster = km$cluster
# Using dplyr
library(dplyr)
correlation_edu_pgs <- pheno %>%
  group_by(km_cluster) %>%
  summarise(correlation = cor(education, pgs_parents, use="complete.obs"))
mean(abs(correlation_edu_pgs$correlation))
#[1] 0.04064645

correlation_bmi_pgs <- pheno %>%
  group_by(km_cluster) %>%
  summarise(correlation = cor(bmi, pgs_parents, use="complete.obs"))
mean(abs(correlation_bmi_pgs$correlation))
#[1] 0.02735008


correlation_townsend_pgs <- pheno %>%
  group_by(km_cluster) %>%
  summarise(correlation = cor(townsend_index, pgs_parents, use="complete.obs"))
mean(abs(correlation_townsend_pgs$correlation))
#[1] 0.02434992


#sensitivity analysis using assessment centres, between cluster correlation
#calculate the means for PGS and towsend index and edu by clusters
means = aggregate(pheno[,c("education","townsend_index","pgs_parents","bmi")], list(pheno$assessment_center), mean, na.rm = TRUE)
cor.test(means$townsend_index,means$pgs_parents)
# Pearson's product-moment correlation
# 
# data:  means$townsend_index and means$pgs_parents
# t = 3.0409, df = 19, p-value = 0.006723
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1865476 0.8050125
# sample estimates:
#       cor 
# 0.5721584 


cor.test(means$education,means$pgs_parents)
#	Pearson's product-moment correlation
# 
# data:  means$education and means$pgs_parents
# t = -6.9703, df = 19, p-value = 1.217e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.9367058 -0.6564439
# sample estimates:
#   cor 
# -0.8478644 

cor.test(means$bmi,means$pgs_parents)
#	Pearson's product-moment correlation
# 
# data:  means$bmi and means$pgs_parents
# t = -5.975, df = 19, p-value = 9.469e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.9190437 -0.5776312
# sample estimates:
#   cor 
# -0.8078707 





#check father alone
pgs_parents = fread("plink/ea_score_parents_pruned.sscore")
load("./family/simulation/UKB/matched_parents_pheno.RData")
load("./family/simulation/UKB/selected_final.RData")
pheno = pheno_father
#pheno=pheno[complete.cases(pheno),]
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

#calculate the means for PGS and towsend index and edu by clusters (between-cluster correlation)
means = aggregate(pheno[,c("education","townsend_index","pgs_parents","bmi")], list(km$cluster), mean, na.rm = TRUE)
cor.test(means$townsend_index,means$pgs_parents)
cor.test(means$bmi,means$pgs_parents)
cor.test(means$education,means$pgs_parents)


#within cluster correlation
pheno$km_cluster = km$cluster
# Using dplyr
library(dplyr)
correlation_edu_pgs <- pheno %>%
  group_by(km_cluster) %>%
  summarise(correlation = cor(education, pgs_parents, use="complete.obs"))
mean(abs(correlation_edu_pgs$correlation))

correlation_bmi_pgs <- pheno %>%
  group_by(km_cluster) %>%
  summarise(correlation = cor(bmi, pgs_parents, use="complete.obs"))
mean(abs(correlation_bmi_pgs$correlation))


correlation_townsend_pgs <- pheno %>%
  group_by(km_cluster) %>%
  summarise(correlation = cor(townsend_index, pgs_parents, use="complete.obs"))
mean(abs(correlation_townsend_pgs$correlation))


#sensitivity analysis using assessment centres, between cluster correlation
#calculate the means for PGS and towsend index and edu by clusters
means = aggregate(pheno[,c("education","townsend_index","pgs_parents","bmi")], list(pheno$assessment_center), mean, na.rm = TRUE)
cor.test(means$townsend_index,means$pgs_parents)
cor.test(means$education,means$pgs_parents)
cor.test(means$bmi,means$pgs_parents)






#check mother alone
pgs_parents = fread("plink/ea_score_parents_pruned.sscore")
load("./family/simulation/UKB/matched_parents_pheno.RData")
load("./family/simulation/UKB/selected_final.RData")
pheno = pheno_mother
#pheno=pheno[complete.cases(pheno),]
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

#calculate the means for PGS and towsend index and edu by clusters (between-cluster correlation)
means = aggregate(pheno[,c("education","townsend_index","pgs_parents","bmi")], list(km$cluster), mean, na.rm = TRUE)
cor.test(means$townsend_index,means$pgs_parents)
cor.test(means$bmi,means$pgs_parents)
cor.test(means$education,means$pgs_parents)


#within cluster correlation
pheno$km_cluster = km$cluster
# Using dplyr
library(dplyr)
correlation_edu_pgs <- pheno %>%
  group_by(km_cluster) %>%
  summarise(correlation = cor(education, pgs_parents, use="complete.obs"))
mean(abs(correlation_edu_pgs$correlation))

correlation_bmi_pgs <- pheno %>%
  group_by(km_cluster) %>%
  summarise(correlation = cor(bmi, pgs_parents, use="complete.obs"))
mean(abs(correlation_bmi_pgs$correlation))


correlation_townsend_pgs <- pheno %>%
  group_by(km_cluster) %>%
  summarise(correlation = cor(townsend_index, pgs_parents, use="complete.obs"))
mean(abs(correlation_townsend_pgs$correlation))


#sensitivity analysis using assessment centres, between cluster correlation
#calculate the means for PGS and towsend index and edu by clusters
means = aggregate(pheno[,c("education","townsend_index","pgs_parents","bmi")], list(pheno$assessment_center), mean, na.rm = TRUE)
cor.test(means$townsend_index,means$pgs_parents)
cor.test(means$education,means$pgs_parents)
cor.test(means$bmi,means$pgs_parents)



#Use mother's BMI and EDU in the simulation study



#Now check children's PGS and mother's BMI to see whether this is any population structure
load("./family/simulation/UKB/results/child_pgs_496.RData")
head(id_parental)
identical(as.character(pheno_mother$ID),as.character(id_parental$X2))
pheno=pheno_mother
pheno$pgs_parents = pgs_c
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

#calculate the means for PGS and towsend index and edu by clusters (between-cluster correlation)
means = aggregate(pheno[,c("education","townsend_index","pgs_parents","bmi")], list(km$cluster), mean, na.rm = TRUE)
cor.test(means$townsend_index,means$pgs_parents)
cor.test(means$bmi,means$pgs_parents)
#	Pearson's product-moment correlation
# 
# data:  means$bmi and means$pgs_parents
# t = -3.4844, df = 98, p-value = 0.0007391
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.4960800 -0.1450543
# sample estimates:
#   cor 
# -0.3320129 
cor.test(means$education,means$pgs_parents)
# Pearson's product-moment correlation
# 
# data:  means$education and means$pgs_parents
# t = -2.4551, df = 98, p-value = 0.01585
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.41739194 -0.04648765
# sample estimates:
#        cor 
# -0.2407079 

#check mean BMI of mother and father
pheno=pheno_mother
pheno$pgs_parents = pgs_c
pheno$bmi_m = selected_final$BMI[match(as.character(pheno$ID),as.character(selected_final$ID))]
pheno$bmi_f = selected_final$BMI[match(as.character(id_parental$X1),as.character(selected_final$ID))]
pheno$bmi = (pheno$bmi_m+pheno$bmi_f)/2
pheno$bmi = scale(pheno$bmi)

tmp = pheno[,2:3]
tmp = scale(tmp)
set.seed(05062024)
km <- kmeans(tmp, 100)
table(km$cluster)

#calculate the means for PGS and towsend index and edu by clusters (between-cluster correlation)
means = aggregate(pheno[,c("pgs_parents","bmi")], list(km$cluster), mean, na.rm = TRUE)
cor.test(means$bmi,means$pgs_parents)
#	Pearson's product-moment correlation
# 
# data:  means$bmi and means$pgs_parents
# t = -3.9078, df = 98, p-value = 0.0001715
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.5256813 -0.1840294
# sample estimates:
#   cor 
# -0.3671753 


#within cluster correlation
pheno$km_cluster = km$cluster
# Using dplyr
library(dplyr)
correlation_edu_pgs <- pheno %>%
  group_by(km_cluster) %>%
  summarise(correlation = cor(education, pgs_parents, use="complete.obs"))
mean(abs(correlation_edu_pgs$correlation))

correlation_bmi_pgs <- pheno %>%
  group_by(km_cluster) %>%
  summarise(correlation = cor(bmi, pgs_parents, use="complete.obs"))
mean(abs(correlation_bmi_pgs$correlation))


correlation_townsend_pgs <- pheno %>%
  group_by(km_cluster) %>%
  summarise(correlation = cor(townsend_index, pgs_parents, use="complete.obs"))
mean(abs(correlation_townsend_pgs$correlation))


#sensitivity analysis using assessment centres, between cluster correlation
#calculate the means for PGS and towsend index and edu by clusters
means = aggregate(pheno[,c("education","townsend_index","pgs_parents","bmi")], list(pheno$assessment_center), mean, na.rm = TRUE)
cor.test(means$townsend_index,means$pgs_parents)
cor.test(means$education,means$pgs_parents)
cor.test(means$bmi,means$pgs_parents)


