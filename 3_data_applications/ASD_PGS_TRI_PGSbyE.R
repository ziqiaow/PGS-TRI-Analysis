#ASD Data Analysis
#PGS-TRI for PGSxE interactions
#April 24, 2024


setwd("./family/data/autism")
source("/users/zwang4/family/R/PGS-TRI.R")

data.model <- readRDS("clean_data_asd_input.rds")


# Fit PRS x E using PRS-TRI
res_int_all=list()
data.model1=data.model[,c(5,65:67,13,14,17,20,23,25,26,29,40,43:46,36,68)]
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

data.model1$sex=data.model1$sex-1
#data.model1=data.model[,c(5,65:67,13,14,17,20,23,25,26,29,40,43:46,36,68)]
#dim(data.model1)
h=0
name=c('asthma','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','sex')
for(i in name){
  h=h+1
  res_int_all[[h]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_combine, 
                           pgs_mother = data.model1$mother_prs_normalized_1000g_combine, 
                           pgs_father = data.model1$father_prs_normalized_1000g_combine,
                           GxE_int = TRUE,
                           formula = as.formula(paste("~",i)),
                           E = data.model1, 
                           side = 2)$res_beta
  
  
}

h=h+1
res_int_all[[h]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_combine, 
                         pgs_mother = data.model1$mother_prs_normalized_1000g_combine, 
                         pgs_father = data.model1$father_prs_normalized_1000g_combine,
                         GxE_int = TRUE,
                         formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                         E = data.model1, 
                         side = 2)$res_beta

res_int=list()

(res_int[[1]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_combine,
                      GxE_int = TRUE,
                      formula = ~ asthma+pre_vitamin+fever+gdiabetes+hyperemesis+pre_labor+eclampsia_preclam+depression_mentalill+child_bwt_low+factor(smoke_tob_ever_preg)+factor(alcohol_ever_preg),
                      E = data.model1, 
                      side = 2)$res_beta)



#We observed significant PRS x Maternal alcohol consumption for two groups: ever 1/preg 0 and ever 1/preg 1 compared to ever 0 (those who never drink alcohol), with a larger OR for those who have ever drink alcohol and who drank alcohol during pregnancy. Next let's check if there is any dose-response effect for alcohol intake frequency during pregnancy.
#Alcohol consumption:
#1.Never
#2.Ever, No during pregnancy
#3.Ever, Yes during pregnancy


res_int_freq=list()
#add alcohol frequency into the groups
data.model1$alcohol_ever1_preg0 <- ifelse(data.model1$alcohol_ever_preg == 1, 1, 0)
data.model1$alcohol_ever1_preg1_mild <- ifelse( (data.model1$alcohol_ever_preg == 2 & data.model1$alcohol_frequency ==4) , 1, 0)
data.model1$alcohol_ever1_preg1_moderate <- ifelse( (data.model1$alcohol_ever_preg == 2 & data.model1$alcohol_frequency ==3) , 1, 0)
data.model1$alcohol_ever1_preg1_heavy <- ifelse( (data.model1$alcohol_ever_preg == 2 & (data.model1$alcohol_frequency ==1 | data.model1$alcohol_frequency ==2)) , 1, 0)
(res_int_freq[[1]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_combine, 
                           pgs_mother = data.model1$mother_prs_normalized_1000g_combine, 
                           pgs_father = data.model1$father_prs_normalized_1000g_combine,
                           GxE_int = TRUE,
                           formula = ~ alcohol_ever1_preg0+alcohol_ever1_preg1_mild+alcohol_ever1_preg1_moderate+alcohol_ever1_preg1_heavy,
                           E = data.model1, 
                           side = 2)$res_beta)

#Higher frequency of maternal alcohol consumption during pregnancy also led to greater interaction effect for mild and moderate frequency drinkers. 


####################################################################
####################################################################
####################################################################

# Remove AMR group in the analysis
data.model2=data.model1[which(data.model1$race != "AMR"),]


res_int_noAMR=list()
table(data.model2$edu_mom)

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
data.model2$edu_mom1=data.model2$edu_mom
data.model2$edu_mom1[which(data.model2$edu_mom == 1 | data.model2$edu_mom == 2 | data.model2$edu_mom == 3 | data.model2$edu_mom == 4 | data.model2$edu_mom == 5 )]=0
data.model2$edu_mom1[which(data.model2$edu_mom == 6 )]=1
data.model2$edu_mom1[which(data.model2$edu_mom == 7 |  data.model2$edu_mom == 8)]=2

#create dummy variables for edu_mom1
data.model2$edu_mom1_bachlor <- ifelse(data.model2$edu_mom1 == 1, 1, 0)
data.model2$edu_mom1_advanced <- ifelse(data.model2$edu_mom1 == 2, 1, 0)

#categorize mother's age at birth
data.model2$mom_age_birth_cate <- cut(data.model2$dob_mom, 
                                      breaks=c(-Inf,20,25,30,35,40,Inf), 
                                      labels=c("<=20","21-25","26-30","31-35","36-40",">40"))

data.model2$sex=data.model2$sex-1
#data.model2=data.model[,c(5,65:67,13,14,17,20,23,25,26,29,40,43:46,36,68)]
#dim(data.model2)
h=0
name=c('asthma','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','sex')
for(i in name){
  h=h+1
  res_int_noAMR[[h]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                             pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                             pgs_father = data.model2$father_prs_normalized_1000g_combine,
                             GxE_int = TRUE,
                             formula = as.formula(paste("~",i)),
                             E = data.model2, 
                             side = 2)$res_beta
  
  
}
h=h+1
res_int_noAMR[[h]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                           pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                           pgs_father = data.model2$father_prs_normalized_1000g_combine,
                           GxE_int = TRUE,
                           formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                           E = data.model2, 
                           side = 2)$res_beta

(res_int[[2]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                      pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                      pgs_father = data.model2$father_prs_normalized_1000g_combine,
                      GxE_int = TRUE,
                      formula = ~ asthma+pre_vitamin+fever+gdiabetes+hyperemesis+pre_labor+eclampsia_preclam+depression_mentalill+child_bwt_low+factor(smoke_tob_ever_preg)+factor(alcohol_ever_preg),
                      E = data.model2, 
                      side = 2)$res_beta)


#add alcohol frequency into the groups
data.model2$alcohol_ever1_preg0 <- ifelse(data.model2$alcohol_ever_preg == 1, 1, 0)
data.model2$alcohol_ever1_preg1_mild <- ifelse( (data.model2$alcohol_ever_preg == 2 & data.model2$alcohol_frequency ==4) , 1, 0)
data.model2$alcohol_ever1_preg1_moderate <- ifelse( (data.model2$alcohol_ever_preg == 2 & data.model2$alcohol_frequency ==3) , 1, 0)
data.model2$alcohol_ever1_preg1_heavy <- ifelse( (data.model2$alcohol_ever_preg == 2 & (data.model2$alcohol_frequency ==1 | data.model2$alcohol_frequency ==2)) , 1, 0)
(res_int_freq[[2]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                           pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                           pgs_father = data.model2$father_prs_normalized_1000g_combine,
                           GxE_int = TRUE,
                           formula = ~ alcohol_ever1_preg0+alcohol_ever1_preg1_mild+alcohol_ever1_preg1_moderate+alcohol_ever1_preg1_heavy,
                           E = data.model2, 
                           side = 2)$res_beta)






#################################################################
#################################################################
#################################################################

# Check the interaction effect in each ancestry group

## European Ancestry
data.model2=data.model1[which(data.model1$race=="EUR"),]
res_int_EUR=list()
h=0
name=c('asthma','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','sex')
for(i in name){
  h=h+1
  res_int_EUR[[h]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                           pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                           pgs_father = data.model2$father_prs_normalized_1000g_combine,
                           GxE_int = TRUE,
                           formula = as.formula(paste("~",i)),
                           E = data.model2, 
                           side = 2)$res_beta
  
  
}
h=h+1
res_int_EUR[[h]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                         pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                         pgs_father = data.model2$father_prs_normalized_1000g_combine,
                         GxE_int = TRUE,
                         formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                         E = data.model2, 
                         side = 2)$res_beta

(res_int[[3]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                      pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                      pgs_father = data.model2$father_prs_normalized_1000g_combine,
                      GxE_int = TRUE,
                      formula = ~ asthma+pre_vitamin+fever+gdiabetes+hyperemesis+pre_labor+eclampsia_preclam+depression_mentalill+child_bwt_low+factor(smoke_tob_ever_preg)+factor(alcohol_ever_preg),
                      E = data.model2, 
                      side = 2)$res_beta)

#--------------------------------------
#Add alcohol frequency

(res_int_freq[[3]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                           pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                           pgs_father = data.model2$father_prs_normalized_1000g_combine,
                           GxE_int = TRUE,
                           formula = ~ alcohol_ever1_preg0+alcohol_ever1_preg1_mild+alcohol_ever1_preg1_moderate+alcohol_ever1_preg1_heavy,
                           E = data.model2, 
                           side = 2)$res_beta)




#################################################################
#################################################################
#################################################################

## African Ancestry
data.model2=data.model1[which(data.model1$race=="AFR"),]

res_int_AFR=list()
h=0
name=c('asthma','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','sex')
for(i in name){
  h=h+1
  res_int_AFR[[h]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                           pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                           pgs_father = data.model2$father_prs_normalized_1000g_combine,
                           GxE_int = TRUE,
                           formula = as.formula(paste("~",i)),
                           E = data.model2, 
                           side = 2,  smalltriosize = T)$res_beta
  
  
}
h=h+1
res_int_AFR[[h]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                         pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                         pgs_father = data.model2$father_prs_normalized_1000g_combine,
                         GxE_int = TRUE,
                         formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                         E = data.model2, 
                         side = 2, smalltriosize = T)$res_beta

(res_int[[4]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                      pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                      pgs_father = data.model2$father_prs_normalized_1000g_combine,
                      GxE_int = TRUE,
                      formula = ~ asthma+pre_vitamin+fever+gdiabetes+hyperemesis+pre_labor+eclampsia_preclam+depression_mentalill+child_bwt_low+factor(smoke_tob_ever_preg)+factor(alcohol_ever_preg),
                      E = data.model2, 
                      side = 2,  smalltriosize = T)$res_beta)

#--------------------------------------
#Add alcohol frequency

(res_int_freq[[4]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                           pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                           pgs_father = data.model2$father_prs_normalized_1000g_combine,
                           GxE_int = TRUE,
                           formula = ~ alcohol_ever1_preg0+alcohol_ever1_preg1_mild+alcohol_ever1_preg1_moderate+alcohol_ever1_preg1_heavy,
                           E = data.model2, 
                           side = 2,  smalltriosize = T)$res_beta)





################################################################
################################################################
################################################################
## Americas
data.model2=data.model1[which(data.model1$race=="AMR"),]
#--------------------------------------
res_int_AMR=list()
h=0
name=c('asthma','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','sex')
for(i in name){
  h=h+1
  res_int_AMR[[h]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                           pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                           pgs_father = data.model2$father_prs_normalized_1000g_combine,
                           GxE_int = TRUE,
                           formula = as.formula(paste("~",i)),
                           E = data.model2, 
                           side = 2)$res_beta
  
  
}
h=h+1
res_int_AMR[[h]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                         pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                         pgs_father = data.model2$father_prs_normalized_1000g_combine,
                         GxE_int = TRUE,
                         formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                         E = data.model2, 
                         side = 2)$res_beta


(res_int[[5]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                      pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                      pgs_father = data.model2$father_prs_normalized_1000g_combine,
                      GxE_int = TRUE,
                      formula = ~ asthma+pre_vitamin+fever+gdiabetes+hyperemesis+pre_labor+eclampsia_preclam+depression_mentalill+child_bwt_low+factor(smoke_tob_ever_preg)+factor(alcohol_ever_preg),
                      E = data.model2, 
                      side = 2)$res_beta)

#--------------------------------------
#Add alcohol frequency

(res_int_freq[[5]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                           pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                           pgs_father = data.model2$father_prs_normalized_1000g_combine,
                           GxE_int = TRUE,
                           formula = ~ alcohol_ever1_preg0+alcohol_ever1_preg1_mild+alcohol_ever1_preg1_moderate+alcohol_ever1_preg1_heavy,
                           E = data.model2, 
                           side = 2)$res_beta)



################################################################
################################################################
################################################################
## Asian Ancestry
#Combine South Asian and East Asian together due to small trio sizes
data.model2=data.model1[which(data.model1$race=="SAS"| data.model1$race == "EAS"),]
#--------------------------------------
res_int_AS=list()
h=0
name=c('asthma','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','factor(smoke_tob_ever_preg)','factor(alcohol_ever_preg)','factor(mom_age_birth_cate)','sex')
for(i in name){
  h=h+1
  res_int_AS[[h]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                          pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                          pgs_father = data.model2$father_prs_normalized_1000g_combine,
                          GxE_int = TRUE,
                          formula = as.formula(paste("~",i)),
                          E = data.model2, 
                          side = 2,  smalltriosize = T)$res_beta
  
  
}
h=h+1
res_int_AS[[h]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                        pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                        pgs_father = data.model2$father_prs_normalized_1000g_combine,
                        GxE_int = TRUE,
                        formula = ~edu_mom1_bachlor+edu_mom1_advanced,
                        E = data.model2, 
                        side = 2,  smalltriosize = T)$res_beta

(res_int[[6]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                      pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                      pgs_father = data.model2$father_prs_normalized_1000g_combine,
                      GxE_int = TRUE,
                      formula = ~ asthma+pre_vitamin+fever+gdiabetes+hyperemesis+pre_labor+eclampsia_preclam+depression_mentalill+child_bwt_low+factor(smoke_tob_ever_preg)+factor(alcohol_ever_preg),
                      E = data.model2, 
                      side = 2, smalltriosize = T)$res_beta)

#--------------------------------------
#Add alcohol frequency

(res_int_freq[[6]]=PGS.TRI(pgs_offspring = data.model2$offspring_prs_normalized_1000g_combine, 
                           pgs_mother = data.model2$mother_prs_normalized_1000g_combine, 
                           pgs_father = data.model2$father_prs_normalized_1000g_combine,
                           GxE_int = TRUE,
                           formula = ~ alcohol_ever1_preg0+alcohol_ever1_preg1_mild+alcohol_ever1_preg1_moderate+alcohol_ever1_preg1_heavy,
                           E = data.model2, 
                           side = 2, smalltriosize = T)$res_beta)





names(res_int_freq)=c("all","all_no_AMR","EUR","AFR","AMR","SAS/EAS")
names(res_int)=c("all","all_no_AMR","EUR","AFR","AMR","SAS/EAS")

save(res_int_freq,res_int_EUR,res_int_AFR,res_int_noAMR,res_int_AS,res_int_AMR,res_int_all,res_int,file="interaction_results.RData")



#Summarize the results
#interactions
tmp_interaction_noamr=res_int_noAMR[[1]][-1,c(1,2,4)]
for(i in 2:length(res_int_noAMR)){
  tmp_interaction_noamr=rbind(tmp_interaction_noamr,res_int_noAMR[[i]][-1,c(1,2,4)])
  
}
rownames(tmp_interaction_noamr)[c(1:9,19)]=c('asthma','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','sex')
tmp_interaction_noamr=rbind(tmp_interaction_noamr,res_int_freq$all_no_AMR[-1,c(1,2,4)])


tmp_interaction_EUR=res_int_EUR[[1]][-1,c(1,2,4)]
for(i in 2:length(res_int_EUR)){
  tmp_interaction_EUR=rbind(tmp_interaction_EUR,res_int_EUR[[i]][-1,c(1,2,4)])
  
}
rownames(tmp_interaction_EUR)[c(1:9,19)]=c('asthma','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','sex')
tmp_interaction_EUR=rbind(tmp_interaction_EUR,res_int_freq$EUR[-1,c(1,2,4)])


tmp_interaction_AFR=res_int_AFR[[1]][-1,c(1,2,4)]
for(i in 2:length(res_int_AFR)){
  tmp_interaction_AFR=rbind(tmp_interaction_AFR,res_int_AFR[[i]][-1,c(1,2,4)])
  
}
rownames(tmp_interaction_AFR)[c(1:9,19)]=c('asthma','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','sex')
tmp_interaction_AFR=rbind(tmp_interaction_AFR,res_int_freq$AFR[-1,c(1,2,4)])



tmp_interaction_AS=res_int_AS[[1]][-1,c(1,2,4)]
for(i in 2:length(res_int_AS)){
  tmp_interaction_AS=rbind(tmp_interaction_AS,res_int_AS[[i]][-1,c(1,2,4)])
  
}
rownames(tmp_interaction_AS)[c(1:9,10,17)]=c('asthma','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','smoke_ever','sex')
tmp_interaction_AS=rbind(tmp_interaction_AS,res_int_freq$`SAS/EAS`[-1,c(1,2,4)])

save(tmp_interaction_noamr,tmp_interaction_EUR,tmp_interaction_AFR,tmp_interaction_AS,file="results/interaction_results_joint.RData")

res_all=data.frame(rbind(tmp_interaction_noamr,tmp_interaction_EUR,tmp_interaction_AFR,tmp_interaction_AS))
res_all$ancestry=c(rep("All_noAMR",dim(tmp_interaction_noamr)[1]),
                   rep("EUR",dim(tmp_interaction_EUR)[1]),
                   rep("AFR",dim(tmp_interaction_AFR)[1]),
                   rep("SAS/EAS",dim(tmp_interaction_AS)[1]))
write.csv(res_all,file="results/interaction_results_joint.csv")




#summarize the separate interaction results
#res_int_EUR,res_int_AFR,res_int_noAMR,res_int_AS,res_int_AMR
tmp_interaction_EUR = do.call(rbind, lapply(res_int_EUR, function(x) x[-1,]))
rownames(tmp_interaction_EUR)[c(1:9,19)] = c('asthma','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','sex')

tmp_interaction_noamr = do.call(rbind,lapply(res_int_noAMR, function(x) x[-1,]))
rownames(tmp_interaction_noamr)[c(1:9,19)] = c('asthma','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','sex')

tmp_interaction_AFR = do.call(rbind, lapply(res_int_AFR, function(x) x[-1,]))
rownames(tmp_interaction_AFR)[c(1:9,19)] = c('asthma','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','sex')

tmp_interaction_AS = do.call(rbind,lapply(res_int_AS, function(x) x[-1,]))
rownames(tmp_interaction_AS)[c(1:10,17)] = c('asthma','pre_vitamin','fever','gdiabetes','hyperemesis','pre_labor','eclampsia_preclam','depression_mentalill','child_bwt_low','smoke_ever','sex')


save(tmp_interaction_noamr,tmp_interaction_EUR,tmp_interaction_AFR,tmp_interaction_AS,file="results/interaction_results_separate.RData")

res_all=data.frame(rbind(tmp_interaction_noamr,tmp_interaction_EUR,tmp_interaction_AFR,tmp_interaction_AS))
res_all$ancestry=c(rep("All_noAMR",dim(tmp_interaction_noamr)[1]),
                   rep("EUR",dim(tmp_interaction_EUR)[1]),
                   rep("AFR",dim(tmp_interaction_AFR)[1]),
                   rep("SAS/EAS",dim(tmp_interaction_AS)[1]))
write.csv(res_all,file="results/interaction_results_separate.csv")
