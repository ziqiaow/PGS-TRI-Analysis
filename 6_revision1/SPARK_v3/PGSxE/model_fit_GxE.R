#Data Application: PRS x Maternal Environmental Exposure Investigation of Autism
#03/05/2025
#Supplementary table S6
#-------------------------------------------------------------------------------
setwd("~/work/family/data/autism/v3/GxE")

source("~/work/family/R/PGS.TRI/R/PGS-TRI.R")
source("~/work/family/R/PGS.TRI/R/pTDT.R")

#load data and data cleaning
library(data.table)
library(fastDummies)
load("~/work/family/data/autism/v3/GxE/autism_model_data.RData")

#combine smoke and tobacco
pheno_clean$smoke_tobacco=pheno_clean$smoke_ever
pheno_clean$smoke_tobacco[which(pheno_clean$tob_ever ==1 )]=1

pheno_clean$smoke_tob_preg=pheno_clean$smoke_preg
pheno_clean$smoke_tob_preg[which(pheno_clean$tob_preg ==1)]=1

#Combine depression and mental illness
pheno_clean$depression_mentalill=pheno_clean$depression
pheno_clean$depression_mentalill[which(pheno_clean$mental_illness==1)]=1

#Categorize the child weight to a binary variable given the small sample size (whether it is low birth weight or not)
pheno_clean$child_bwt_low=pheno_clean$child_bwt_g
pheno_clean$child_bwt_low[which(pheno_clean$child_bwt_g >= 2500)] = 0
pheno_clean$child_bwt_low[which(pheno_clean$child_bwt_g < 2500)] = 1

#create a new variable: alcohol: ever/preg: 0/0,1/0,1/1
#smoke_tobacco: ever/preg: 0/0,1/0,1/1
pheno_clean$alcohol_ever_preg=pheno_clean$alcohol_ever #group 0: 0/0
pheno_clean$alcohol_ever_preg[which(pheno_clean$alcohol_pregnancy == 0 & pheno_clean$alcohol_ever == 1)] = 1
pheno_clean$alcohol_ever_preg[which(pheno_clean$alcohol_pregnancy == 1 & pheno_clean$alcohol_ever == 1)] = 2

pheno_clean$smoke_tob_ever_preg=pheno_clean$smoke_tobacco #group 0: 0/0
pheno_clean$smoke_tob_ever_preg[which(pheno_clean$smoke_tob_preg == 0 & pheno_clean$smoke_tobacco == 1)] = 1
pheno_clean$smoke_tob_ever_preg[which(pheno_clean$smoke_tob_preg == 1 & pheno_clean$smoke_tobacco == 1)] = 2

#add previously standardized PGS values to trio phenotype data
#raw PGS
pheno_clean$pgs_child = PGS_pheno$SCORE1_SUM[match(pheno_clean$Subject_ParticipantId,PGS_pheno$IID)]
pheno_clean$pgs_mother = PGS_pheno$SCORE1_SUM[match(pheno_clean$Authorizer_ParticipantId,PGS_pheno$IID)]
pheno_clean$pgs_father = PGS_pheno$SCORE1_SUM[match(pheno_clean$biofather_sp_id,PGS_pheno$IID)]

#PC PGS
pheno_clean$pgs_child_PC = PGS_pheno$ASD_PGS_adj_mean[match(pheno_clean$Subject_ParticipantId,PGS_pheno$IID)]
pheno_clean$pgs_mother_PC = PGS_pheno$ASD_PGS_adj_mean[match(pheno_clean$Authorizer_ParticipantId,PGS_pheno$IID)]
pheno_clean$pgs_father_PC = PGS_pheno$ASD_PGS_adj_mean[match(pheno_clean$biofather_sp_id,PGS_pheno$IID)]

#1KG mean/sd PGS
pheno_clean$pgs_child_1KG = PGS_pheno$ASD_PGS_1KGstand[match(pheno_clean$Subject_ParticipantId,PGS_pheno$IID)]
pheno_clean$pgs_mother_1KG = PGS_pheno$ASD_PGS_1KGstand[match(pheno_clean$Authorizer_ParticipantId,PGS_pheno$IID)]
pheno_clean$pgs_father_1KG = PGS_pheno$ASD_PGS_1KGstand[match(pheno_clean$biofather_sp_id,PGS_pheno$IID)]






#Check homogeneous family
pheno_homo = pheno_clean[which(pheno_clean$ancestry == pheno_clean$ancestry_f & pheno_clean$ancestry == pheno_clean$ancestry_m),]
pheno_homo = pheno_homo[which(pheno_homo$ancestry %in% c("AFR", "AMR", "EUR", "EAS", "SAS")),]
#1748 families
table(pheno_homo$ancestry)
# AFR      AMR      EAS      EUR      SAS 
#51       92       37     1543       25       

# Fit PGS x E using PGS-TRI

data.model1 = pheno_homo
table(data.model1$edu_mom)

#1= Less than high school
# 2 = 12 years, completed high school or equivalent
# 3 = 1 – 3 years of college
# 4 = Completed technical college
# 5 = Associate’s degree
# 6 = 4 years of college or bachelor’s degree
# 7 = Master’s degree
# 8 = Advanced degree
# 888 = Don’t know
# 999 = Prefer not to answer

#Group anyone who does not have bachelor's degree to one group, master's and advanced degree to one group
data.model1$edu_mom1=data.model1$edu_mom
data.model1$edu_mom1[which(data.model1$edu_mom == 1 | data.model1$edu_mom == 2 | data.model1$edu_mom == 3 | data.model1$edu_mom == 4 | data.model1$edu_mom == 5 )]=0
data.model1$edu_mom1[which(data.model1$edu_mom == 6 )]=1
data.model1$edu_mom1[which(data.model1$edu_mom == 7 |  data.model1$edu_mom == 8)]=2

#create dummy variables for edu_mom1
data.model1$edu_mom1_bachlor <- ifelse(data.model1$edu_mom1 == 1, 1, 0)
data.model1$edu_mom1_advanced <- ifelse(data.model1$edu_mom1 == 2, 1, 0)

#categorize mother's age at birth
data.model1$mom_age_birth_cate <- cut(data.model1$dob_mom, 
                                      breaks=c(-Inf,20,25,30,35,40,Inf), 
                                      labels=c("<=20","21-25","26-30","31-35","36-40",">40"))



# #add alcohol frequency into the groups
# data.model1$alcohol_ever1_preg0 <- ifelse(data.model1$alcohol_ever_preg == 1, 1, 0)
# data.model1$alcohol_ever1_preg1_mild <- ifelse( (data.model1$alcohol_ever_preg == 2 & data.model1$alcohol_frequency ==4) , 1, 0)
# data.model1$alcohol_ever1_preg1_moderate <- ifelse( (data.model1$alcohol_ever_preg == 2 & data.model1$alcohol_frequency ==3) , 1, 0)
# data.model1$alcohol_ever1_preg1_heavy <- ifelse( (data.model1$alcohol_ever_preg == 2 & (data.model1$alcohol_frequency ==1 | data.model1$alcohol_frequency ==2)) , 1, 0)


#All
res_int_all=list()
h=0
name=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','factor(sex)')
for(i in name){
  h=h+1
  res_int_all[[h]]=PGS.TRI(pgs_offspring = data.model1$pgs_child_PC, 
                           pgs_mother = data.model1$pgs_mother_PC, 
                           pgs_father = data.model1$pgs_father_PC,
                           GxE_int = TRUE,
                           formula = as.formula(paste("~",i)),
                           E = data.model1, 
                           side = 2)$res_beta
  
  
}

h=h+1
res_int_all[[h]]=PGS.TRI(pgs_offspring = data.model1$pgs_child_PC, 
                         pgs_mother = data.model1$pgs_mother_PC, 
                         pgs_father = data.model1$pgs_father_PC,
                         GxE_int = TRUE,
                         formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                         E = data.model1, 
                         side = 2)$res_beta

# #--------------------------------------
# #Add alcohol frequency
# h=h+1
# (res_int_all[[h]]=PGS.TRI(pgs_offspring = data.model1$pgs_child_PC, 
#                            pgs_mother = data.model1$pgs_mother_PC, 
#                            pgs_father = data.model1$pgs_father_PC,
#                            GxE_int = TRUE,
#                            formula = ~ alcohol_ever1_preg0+alcohol_ever1_preg1_mild+alcohol_ever1_preg1_moderate+alcohol_ever1_preg1_heavy,
#                            E = data.model1, 
#                            side = 2)$res_beta)



h=h+1
(res_int_all[[h]]=PGS.TRI(pgs_offspring = data.model1$pgs_child_PC, 
                      pgs_mother = data.model1$pgs_mother_PC, 
                      pgs_father = data.model1$pgs_father_PC,
                      GxE_int = TRUE,
                      formula = ~ asthma+pre_vitamin+fever+gdiabetes+hyperemesis+pre_labor+eclampsia_preclam+depression_mentalill+child_bwt_low+factor(smoke_tob_ever_preg)+factor(alcohol_ever_preg),
                      E = data.model1, 
                      side = 2)$res_beta)




#EUR
data.model2 = data.model1[which(data.model1$ancestry == "EUR"),]
res_int_EUR=list()
h=0
name=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','factor(sex)')
for(i in name){
  h=h+1
  res_int_EUR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                           pgs_mother = data.model2$pgs_mother_PC, 
                           pgs_father = data.model2$pgs_father_PC,
                           GxE_int = TRUE,
                           formula = as.formula(paste("~",i)),
                           E = data.model2, 
                           side = 2)$res_beta
  
  
}

h=h+1
res_int_EUR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                         pgs_mother = data.model2$pgs_mother_PC, 
                         pgs_father = data.model2$pgs_father_PC,
                         GxE_int = TRUE,
                         formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                         E = data.model2, 
                         side = 2)$res_beta

h=h+1
(res_int_EUR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                      pgs_mother = data.model2$pgs_mother_PC, 
                      pgs_father = data.model2$pgs_father_PC,
                      GxE_int = TRUE,
                      formula = ~ asthma+pre_vitamin+fever+gdiabetes+hyperemesis+pre_labor+eclampsia_preclam+depression_mentalill+child_bwt_low+factor(smoke_tob_ever_preg)+factor(alcohol_ever_preg),
                      E = data.model2, 
                      side = 2)$res_beta)





#AFR
data.model2 = data.model1[which(data.model1$ancestry == "AFR"),]
res_int_AFR=list()
h=0
name=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','factor(sex)')
for(i in name){
  h=h+1
  res_int_AFR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                           pgs_mother = data.model2$pgs_mother_PC, 
                           pgs_father = data.model2$pgs_father_PC,
                           GxE_int = TRUE,
                           formula = as.formula(paste("~",i)),
                           E = data.model2, 
                           side = 2)$res_beta
  
  
}

h=h+1
res_int_AFR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                         pgs_mother = data.model2$pgs_mother_PC, 
                         pgs_father = data.model2$pgs_father_PC,
                         GxE_int = TRUE,
                         formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                         E = data.model2, 
                         side = 2)$res_beta



#AMR
data.model2 = data.model1[which(data.model1$ancestry == "AMR"),]
res_int_AMR=list()
h=0
name=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','factor(sex)')
for(i in name){
  h=h+1
  res_int_AMR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                           pgs_mother = data.model2$pgs_mother_PC, 
                           pgs_father = data.model2$pgs_father_PC,
                           GxE_int = TRUE,
                           formula = as.formula(paste("~",i)),
                           E = data.model2, 
                           side = 2)$res_beta
  
  
}

h=h+1
res_int_AMR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                         pgs_mother = data.model2$pgs_mother_PC, 
                         pgs_father = data.model2$pgs_father_PC,
                         GxE_int = TRUE,
                         formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                         E = data.model2, 
                         side = 2)$res_beta

h=h+1
(res_int_AMR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                          pgs_mother = data.model2$pgs_mother_PC, 
                          pgs_father = data.model2$pgs_father_PC,
                          GxE_int = TRUE,
                          formula = ~ asthma+pre_vitamin+fever+gdiabetes+hyperemesis+pre_labor+eclampsia_preclam+depression_mentalill+child_bwt_low+factor(smoke_tob_ever_preg)+factor(alcohol_ever_preg),
                          E = data.model2, 
                          side = 2)$res_beta)

#SAS/EAS
data.model2 = data.model1[which(data.model1$ancestry == "SAS" | data.model1$ancestry == "EAS"),]
res_int_SAS=list()
h=0
name=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','factor(sex)')
for(i in name){
  h=h+1
  res_int_SAS[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                           pgs_mother = data.model2$pgs_mother_PC, 
                           pgs_father = data.model2$pgs_father_PC,
                           GxE_int = TRUE,
                           formula = as.formula(paste("~",i)),
                           E = data.model2, 
                           side = 2)$res_beta
  
  
}

h=h+1
res_int_SAS[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                         pgs_mother = data.model2$pgs_mother_PC, 
                         pgs_father = data.model2$pgs_father_PC,
                         GxE_int = TRUE,
                         formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                         E = data.model2, 
                         side = 2)$res_beta


h=h+1
(res_int_SAS[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                          pgs_mother = data.model2$pgs_mother_PC, 
                          pgs_father = data.model2$pgs_father_PC,
                          GxE_int = TRUE,
                          formula = ~ asthma+pre_vitamin+fever+gdiabetes+hyperemesis+pre_labor+eclampsia_preclam+depression_mentalill+child_bwt_low+factor(smoke_tob_ever_preg)+factor(alcohol_ever_preg),
                          E = data.model2, 
                          side = 2)$res_beta)


save(res_int_EUR,res_int_AFR,res_int_AMR,res_int_all,res_int_SAS,file="homo/interaction_results_v3.RData")

#summarize the results

#interactions
tmp_interaction_all=res_int_all[[1]][-1,c(1,2,4)]
for(i in 2:18){
  tmp_interaction_all=rbind(tmp_interaction_all,res_int_all[[i]][-1,c(1,2,4)])
  
}
rownames(tmp_interaction_all)[c(1:13,23)]=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','sex')



tmp_interaction_EUR=res_int_EUR[[1]][-1,c(1,2,4)]
for(i in 2:18){
  tmp_interaction_EUR=rbind(tmp_interaction_EUR,res_int_EUR[[i]][-1,c(1,2,4)])
  
}
rownames(tmp_interaction_EUR)[c(1:13,23)]=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','sex')


tmp_interaction_AFR=res_int_AFR[[1]][-1,c(1,2,4)]
for(i in 2:18){
  tmp_interaction_AFR=rbind(tmp_interaction_AFR,res_int_AFR[[i]][-1,c(1,2,4)])
  
}
rownames(tmp_interaction_AFR)[c(1:13,22)]=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','sex')

tmp_interaction_AMR=res_int_AMR[[1]][-1,c(1,2,4)]
for(i in 2:18){
  tmp_interaction_AMR=rbind(tmp_interaction_AMR,res_int_AMR[[i]][-1,c(1,2,4)])
  
}
rownames(tmp_interaction_AMR)[c(1:13,23)]=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','sex')


tmp_interaction_SAS=res_int_SAS[[1]][-1,c(1,2,4)]
for(i in 2:18){
  tmp_interaction_SAS=rbind(tmp_interaction_SAS,res_int_SAS[[i]][-1,c(1,2,4)])
  
}
rownames(tmp_interaction_SAS)[c(1:14,20)]=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','smoke','sex')

save(tmp_interaction_all,tmp_interaction_EUR,tmp_interaction_AFR,tmp_interaction_AMR,tmp_interaction_SAS,file="homo/interaction_results_clean.RData")

res_all=data.frame(rbind(tmp_interaction_all,tmp_interaction_EUR,tmp_interaction_AFR,tmp_interaction_AMR,tmp_interaction_SAS))
res_all$ancestry=c(rep("All",dim(tmp_interaction_all)[1]),
                   rep("EUR",dim(tmp_interaction_EUR)[1]),
                   rep("AFR",dim(tmp_interaction_AFR)[1]),
                   rep("AMR",dim(tmp_interaction_AMR)[1]),
                   rep("SAS/EAS",dim(tmp_interaction_SAS)[1]))
write.csv(res_all,file="homo/interaction_results_clean_v3.csv")

res_all[,c(1,2,3)] = round(res_all[,c(1,2,3)],digits = 3)
write.csv(res_all,file="homo/interaction_results_clean_round_v3.csv")





#Check homogeneous+heterogeneous family
pheno_homoheter = pheno_clean[which(pheno_clean$ancestry %in% c("AFR", "AMR", "EUR", "EAS", "SAS")),]
#1748 families
table(pheno_homoheter$ancestry)
# AFR      AMR      EAS      EUR      SAS 
#92  207   40 1635   37 

# Fit PGS x E using PGS-TRI

data.model1 = pheno_homoheter
table(data.model1$edu_mom)

#1= Less than high school
# 2 = 12 years, completed high school or equivalent
# 3 = 1 – 3 years of college
# 4 = Completed technical college
# 5 = Associate’s degree
# 6 = 4 years of college or bachelor’s degree
# 7 = Master’s degree
# 8 = Advanced degree
# 888 = Don’t know
# 999 = Prefer not to answer

#Group anyone who does not have bachelor's degree to one group, master's and advanced degree to one group
data.model1$edu_mom1=data.model1$edu_mom
data.model1$edu_mom1[which(data.model1$edu_mom == 1 | data.model1$edu_mom == 2 | data.model1$edu_mom == 3 | data.model1$edu_mom == 4 | data.model1$edu_mom == 5 )]=0
data.model1$edu_mom1[which(data.model1$edu_mom == 6 )]=1
data.model1$edu_mom1[which(data.model1$edu_mom == 7 |  data.model1$edu_mom == 8)]=2

#create dummy variables for edu_mom1
data.model1$edu_mom1_bachlor <- ifelse(data.model1$edu_mom1 == 1, 1, 0)
data.model1$edu_mom1_advanced <- ifelse(data.model1$edu_mom1 == 2, 1, 0)

#categorize mother's age at birth
data.model1$mom_age_birth_cate <- cut(data.model1$dob_mom, 
                                      breaks=c(-Inf,20,25,30,35,40,Inf), 
                                      labels=c("<=20","21-25","26-30","31-35","36-40",">40"))



# #add alcohol frequency into the groups
# data.model1$alcohol_ever1_preg0 <- ifelse(data.model1$alcohol_ever_preg == 1, 1, 0)
# data.model1$alcohol_ever1_preg1_mild <- ifelse( (data.model1$alcohol_ever_preg == 2 & data.model1$alcohol_frequency ==4) , 1, 0)
# data.model1$alcohol_ever1_preg1_moderate <- ifelse( (data.model1$alcohol_ever_preg == 2 & data.model1$alcohol_frequency ==3) , 1, 0)
# data.model1$alcohol_ever1_preg1_heavy <- ifelse( (data.model1$alcohol_ever_preg == 2 & (data.model1$alcohol_frequency ==1 | data.model1$alcohol_frequency ==2)) , 1, 0)


#All
res_int_all=list()
h=0
name=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','factor(sex)')
for(i in name){
  h=h+1
  res_int_all[[h]]=PGS.TRI(pgs_offspring = data.model1$pgs_child_PC, 
                           pgs_mother = data.model1$pgs_mother_PC, 
                           pgs_father = data.model1$pgs_father_PC,
                           GxE_int = TRUE,
                           formula = as.formula(paste("~",i)),
                           E = data.model1, 
                           side = 2)$res_beta
  
  
}

h=h+1
res_int_all[[h]]=PGS.TRI(pgs_offspring = data.model1$pgs_child_PC, 
                         pgs_mother = data.model1$pgs_mother_PC, 
                         pgs_father = data.model1$pgs_father_PC,
                         GxE_int = TRUE,
                         formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                         E = data.model1, 
                         side = 2)$res_beta

# #--------------------------------------
# #Add alcohol frequency
# h=h+1
# (res_int_all[[h]]=PGS.TRI(pgs_offspring = data.model1$pgs_child_PC, 
#                            pgs_mother = data.model1$pgs_mother_PC, 
#                            pgs_father = data.model1$pgs_father_PC,
#                            GxE_int = TRUE,
#                            formula = ~ alcohol_ever1_preg0+alcohol_ever1_preg1_mild+alcohol_ever1_preg1_moderate+alcohol_ever1_preg1_heavy,
#                            E = data.model1, 
#                            side = 2)$res_beta)



h=h+1
(res_int_all[[h]]=PGS.TRI(pgs_offspring = data.model1$pgs_child_PC, 
                          pgs_mother = data.model1$pgs_mother_PC, 
                          pgs_father = data.model1$pgs_father_PC,
                          GxE_int = TRUE,
                          formula = ~ asthma+pre_vitamin+fever+gdiabetes+hyperemesis+pre_labor+eclampsia_preclam+depression_mentalill+child_bwt_low+factor(smoke_tob_ever_preg)+factor(alcohol_ever_preg),
                          E = data.model1, 
                          side = 2)$res_beta)




#EUR
data.model2 = data.model1[which(data.model1$ancestry == "EUR"),]
res_int_EUR=list()
h=0
name=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','factor(sex)')
for(i in name){
  h=h+1
  res_int_EUR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                           pgs_mother = data.model2$pgs_mother_PC, 
                           pgs_father = data.model2$pgs_father_PC,
                           GxE_int = TRUE,
                           formula = as.formula(paste("~",i)),
                           E = data.model2, 
                           side = 2)$res_beta
  
  
}

h=h+1
res_int_EUR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                         pgs_mother = data.model2$pgs_mother_PC, 
                         pgs_father = data.model2$pgs_father_PC,
                         GxE_int = TRUE,
                         formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                         E = data.model2, 
                         side = 2)$res_beta

h=h+1
(res_int_EUR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                          pgs_mother = data.model2$pgs_mother_PC, 
                          pgs_father = data.model2$pgs_father_PC,
                          GxE_int = TRUE,
                          formula = ~ asthma+pre_vitamin+fever+gdiabetes+hyperemesis+pre_labor+eclampsia_preclam+depression_mentalill+child_bwt_low+factor(smoke_tob_ever_preg)+factor(alcohol_ever_preg),
                          E = data.model2, 
                          side = 2)$res_beta)





#AFR
data.model2 = data.model1[which(data.model1$ancestry == "AFR"),]
res_int_AFR=list()
h=0
name=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','factor(sex)')
for(i in name){
  h=h+1
  res_int_AFR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                           pgs_mother = data.model2$pgs_mother_PC, 
                           pgs_father = data.model2$pgs_father_PC,
                           GxE_int = TRUE,
                           formula = as.formula(paste("~",i)),
                           E = data.model2, 
                           side = 2)$res_beta
  
  
}

h=h+1
res_int_AFR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                         pgs_mother = data.model2$pgs_mother_PC, 
                         pgs_father = data.model2$pgs_father_PC,
                         GxE_int = TRUE,
                         formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                         E = data.model2, 
                         side = 2)$res_beta



#AMR
data.model2 = data.model1[which(data.model1$ancestry == "AMR"),]
res_int_AMR=list()
h=0
name=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','factor(sex)')
for(i in name){
  h=h+1
  res_int_AMR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                           pgs_mother = data.model2$pgs_mother_PC, 
                           pgs_father = data.model2$pgs_father_PC,
                           GxE_int = TRUE,
                           formula = as.formula(paste("~",i)),
                           E = data.model2, 
                           side = 2)$res_beta
  
  
}

h=h+1
res_int_AMR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                         pgs_mother = data.model2$pgs_mother_PC, 
                         pgs_father = data.model2$pgs_father_PC,
                         GxE_int = TRUE,
                         formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                         E = data.model2, 
                         side = 2)$res_beta

h=h+1
(res_int_AMR[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                          pgs_mother = data.model2$pgs_mother_PC, 
                          pgs_father = data.model2$pgs_father_PC,
                          GxE_int = TRUE,
                          formula = ~ asthma+pre_vitamin+fever+gdiabetes+hyperemesis+pre_labor+eclampsia_preclam+depression_mentalill+child_bwt_low+factor(smoke_tob_ever_preg)+factor(alcohol_ever_preg),
                          E = data.model2, 
                          side = 2)$res_beta)

#SAS/EAS
data.model2 = data.model1[which(data.model1$ancestry == "SAS" | data.model1$ancestry == "EAS"),]
res_int_SAS=list()
h=0
name=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','factor(sex)')
for(i in name){
  h=h+1
  res_int_SAS[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                           pgs_mother = data.model2$pgs_mother_PC, 
                           pgs_father = data.model2$pgs_father_PC,
                           GxE_int = TRUE,
                           formula = as.formula(paste("~",i)),
                           E = data.model2, 
                           side = 2)$res_beta
  
  
}

h=h+1
res_int_SAS[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                         pgs_mother = data.model2$pgs_mother_PC, 
                         pgs_father = data.model2$pgs_father_PC,
                         GxE_int = TRUE,
                         formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                         E = data.model2, 
                         side = 2)$res_beta


h=h+1
(res_int_SAS[[h]]=PGS.TRI(pgs_offspring = data.model2$pgs_child_PC, 
                          pgs_mother = data.model2$pgs_mother_PC, 
                          pgs_father = data.model2$pgs_father_PC,
                          GxE_int = TRUE,
                          formula = ~ asthma+pre_vitamin+fever+gdiabetes+hyperemesis+pre_labor+eclampsia_preclam+depression_mentalill+child_bwt_low+factor(smoke_tob_ever_preg)+factor(alcohol_ever_preg),
                          E = data.model2, 
                          side = 2)$res_beta)


save(res_int_EUR,res_int_AFR,res_int_AMR,res_int_all,res_int_SAS,file="homo+heter/interaction_results_v3.RData")

#summarize the results

#interactions
tmp_interaction_all=res_int_all[[1]][-1,c(1,2,4)]
for(i in 2:18){
  tmp_interaction_all=rbind(tmp_interaction_all,res_int_all[[i]][-1,c(1,2,4)])
  
}
rownames(tmp_interaction_all)[c(1:13,23)]=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','sex')



tmp_interaction_EUR=res_int_EUR[[1]][-1,c(1,2,4)]
for(i in 2:18){
  tmp_interaction_EUR=rbind(tmp_interaction_EUR,res_int_EUR[[i]][-1,c(1,2,4)])
  
}
rownames(tmp_interaction_EUR)[c(1:13,23)]=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','sex')


tmp_interaction_AFR=res_int_AFR[[1]][-1,c(1,2,4)]
for(i in 2:18){
  tmp_interaction_AFR=rbind(tmp_interaction_AFR,res_int_AFR[[i]][-1,c(1,2,4)])
  
}
rownames(tmp_interaction_AFR)[c(1:13,22)]=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','sex')

tmp_interaction_AMR=res_int_AMR[[1]][-1,c(1,2,4)]
for(i in 2:18){
  tmp_interaction_AMR=rbind(tmp_interaction_AMR,res_int_AMR[[i]][-1,c(1,2,4)])
  
}
rownames(tmp_interaction_AMR)[c(1:13,23)]=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','sex')


tmp_interaction_SAS=res_int_SAS[[1]][-1,c(1,2,4)]
for(i in 2:18){
  tmp_interaction_SAS=rbind(tmp_interaction_SAS,res_int_SAS[[i]][-1,c(1,2,4)])
  
}
rownames(tmp_interaction_SAS)[c(1:14,20)]=c('asthma','mom_bmi','age_dad_birth','dob_mom','mom_wt_gain_kg','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','smoke','sex')

save(tmp_interaction_all,tmp_interaction_EUR,tmp_interaction_AFR,tmp_interaction_AMR,tmp_interaction_SAS,file="homo+heter/interaction_results_clean.RData")

res_all=data.frame(rbind(tmp_interaction_all,tmp_interaction_EUR,tmp_interaction_AFR,tmp_interaction_AMR,tmp_interaction_SAS))
res_all$ancestry=c(rep("All",dim(tmp_interaction_all)[1]),
                   rep("EUR",dim(tmp_interaction_EUR)[1]),
                   rep("AFR",dim(tmp_interaction_AFR)[1]),
                   rep("AMR",dim(tmp_interaction_AMR)[1]),
                   rep("SAS/EAS",dim(tmp_interaction_SAS)[1]))
write.csv(res_all,file="homo+heter/interaction_results_clean_v3.csv")

res_all[,c(1,2,3)] = round(res_all[,c(1,2,3)],digits = 3)
write.csv(res_all,file="homo+heter/interaction_results_clean_round_v3.csv")


