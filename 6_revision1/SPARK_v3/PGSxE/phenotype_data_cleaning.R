#----------------------------------
#phenotype data cleaning
#GEARs v3
#March 5, 2025
#Environmental variables v3:
#/dcs05/ladd/NDEpi/data/MasterCohortData/GEARS/SPARK/iWES_v3/RM0180Fallin_GEARS2/Mothers/GEARS2 Survey_Responses_v5.xlsx
#Environmental variables v1:
#/dcs05/ladd/NDEpi/data/MasterCohortData/GEARS/SPARK/RM0180Fallin_GEARS2/Aug_2022_correct/Mothers/GEARS2 Survey_Responses_v4.xlsx
#1.	Raw genetic data from SPARK live in this directory: /dcs05/ladd/NDEpi/data/MasterCohortData/GEARS/SPARK/iWES_v1 and I attached the README file and the data release notes (also available in the directory). The cleaned & imputed genotype data (only the subset of family members participating in GEARS and with geno data available were imputed, N= 7,121) are also available.
#2.	The GEARS survey (data dictionary attached) was completed in two waves and the raw data for a total of 2,748 mother-child pairs (plus dads for almost all of these) live in two sub-directories within /dcs05/ladd/NDEpi/data/MasterCohortData/GEARS/SPARK: /RM0002Fallin_GEARS and /RM0180Fallin_GEARS2. These survey data (N = 2,748 probands) were hamonized into one dataset, which is available here: /dcs05/ladd/NDEpi/data/Projects/InProgress/lgrosven2/Dissertation/general/data/gears_surv/alldata_July22.csv

#----------------------------------
#https://choishingwan.github.io/PRS-Tutorial/base/
library(data.table)
setwd("~/work/family/data/autism/v3/GxE")
load("~/work/family/data/autism/v3/PRS_SPARK_PCstandardized.RData")
load("~/work/family/data/autism/v3/roles_id_complete_casetrio_ancestry.RData")
#load new phenotype data v3
pheno_v3 = fread("pheno_raw.csv")
any(duplicated(pheno_v3$Authorizer_ParticipantId)) #FALSE
#phenotype data from batch 1
pheno_v1=read.csv("../../documents/alldata_July22.csv",stringsAsFactors = F)
identical(colnames(pheno_v3),colnames(pheno_v1))
#[1] TRUE
any(duplicated(pheno_v1$Authorizer_ParticipantId)) #TRUE
which(duplicated(pheno_v1$Authorizer_ParticipantId)) #this is because each mother have more than 1 children, so different age at birth is recorded for each child, same mom
#[1] 1539 1549 1679 1848 1872 1887 2007 2024 2094 2191 2270 2468 2558 2734
pheno_v1_merge = merge(pheno_v1, id_final, by.x = c("Authorizer_ParticipantId","Subject_ParticipantId"),by.y = c("biomother_sp_id","subject_sp_id"))
any(duplicated(pheno_v1_merge$Authorizer_ParticipantId)) #[1] FALSE
any(duplicated(pheno_v1_merge$Subject_ParticipantId)) #[1] FALSE

pheno_v3_merge = merge(pheno_v3, id_final, by.x = c("Authorizer_ParticipantId","Subject_ParticipantId"),by.y = c("biomother_sp_id","subject_sp_id"))
any(duplicated(pheno_v3_merge$Authorizer_ParticipantId)) #[1] FALSE
any(duplicated(pheno_v3_merge$Subject_ParticipantId)) #[1] FALSE
id_c = intersect(pheno_v1_merge$Subject_ParticipantId,pheno_v3_merge$Subject_ParticipantId)
id_m = intersect(pheno_v1_merge$Authorizer_ParticipantId,pheno_v3_merge$Authorizer_ParticipantId)

#639 duplicated mothers between v1 and v3
#let's use the newest data in v3 for those duplicated moms
#note that for Authorizer_AgeAtSubmission and Subject_AgeAtSubmission v1 is measured by month, v3 measured by year
pheno_v1_merge = pheno_v1_merge[-match(id_c,pheno_v1_merge$Subject_ParticipantId),]
intersect(pheno_v1_merge$Subject_ParticipantId,pheno_v3_merge$Subject_ParticipantId)
intersect(pheno_v1_merge$Authorizer_ParticipantId,pheno_v3_merge$Authorizer_ParticipantId)

#combine both v1 and v3 data
pheno_clean = rbind(pheno_v3_merge,pheno_v1_merge) #2089 mothers
any(duplicated(pheno_clean$Subject_ParticipantId)) #[1] FALSE

id = c(pheno_clean$biofather_sp_id,pheno_clean$Authorizer_ParticipantId,pheno_clean$Subject_ParticipantId)
PGS_pheno = dat[match(id,dat$IID),]
save(pheno_clean,PGS_pheno,file="pheno_PGS_clean.RData")





#----------------------------------------------------------------------------------
#Clean the data for downstream analysis
#do a chi-sq test for depression and mental illness
psychiatric=pheno_clean[,c("depression","mental_illness")]
psychiatric=psychiatric[which(psychiatric$depression==0 | psychiatric$depression==1),]
psychiatric=psychiatric[which(psychiatric$mental_illness!= 999),]
chisq.test(psychiatric$depression,psychiatric$mental_illness)
#	Pearson's Chi-squared test

#data:  psychiatric$depression and psychiatric$mental_illness
#X-squared = 136.05, df = 2, p-value < 2.2e-16

#they are highly corrected, I would suggest to combine these two (test both, separately and combined)

#select: mother's edu, mother's age at birth, pre_vitamin, asthma,pre_labor, depression, mental illness,flu, fever,elampsia,gdiabetes,hyperemesis,preclam,smoke_ever,tob_ever (combine with smoke), smoke_preg/tob_preg, alcohol_preg, alcohol_freq
#child's weight at birth, mother's weight, mother's height
#convert weight from lb to kg
pheno_clean$mom_wt_1[which(pheno_clean$mom_wt_1 == -666 | pheno_clean$mom_wt_1 == 888)]=NA
pheno_clean$mom_wt_2[which(pheno_clean$mom_wt_2 == -666 | pheno_clean$mom_wt_2 == 888)]=NA
pheno_clean$mom_wt_kg=pheno_clean$mom_wt_2
pheno_clean$mom_wt_kg[which(pheno_clean$mom_wt_select == 1)] = pheno_clean$mom_wt_1[which(pheno_clean$mom_wt_select == 1)]* 0.45359237

summary(pheno_clean$mom_wt_kg)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  32.66   58.97   70.00   74.92   85.28  185.52     197 
#convert height from ft/inches to m
#https://www.albireo.ch/bodyconverter/formula.html
#cm=(feet * 30.48) + (inches * 2.54)

pheno_clean$ht_select_2[which(pheno_clean$ht_select_2 == -666 | pheno_clean$ht_select_2 == 888)]=NA
pheno_clean$ht_select_1_ft[which(pheno_clean$ht_select_1_ft == -666 | pheno_clean$ht_select_1_ft == 888)]=NA
pheno_clean$ht_select_1_in[which(pheno_clean$ht_select_1_in == -666 | pheno_clean$ht_select_1_in == 888 | pheno_clean$ht_select_1_in == -999)]=NA
pheno_clean$mom_ht_cm=pheno_clean$ht_select_2
tmp=pheno_clean$ht_select_1_ft * 30.48 + pheno_clean$ht_select_1_in * 2.54
pheno_clean$mom_ht_cm[which(pheno_clean$ht_select == 1)] = tmp[which(pheno_clean$ht_select == 1)]


summary(pheno_clean$mom_ht_cm)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  66.0   160.0   165.1   165.0   170.2   193.0      81 


#calculate mom's BMI
#BMI = weight (kg) / [height (m)]2
pheno_clean$mom_bmi=pheno_clean$mom_wt_kg / (pheno_clean$mom_ht_cm/100)^2
#there is an outlier where the mother input their height in cm to be 66 (might be a typo that missed the 166 or the mother intended to input 66inch). Need to drop this entry.
pheno_clean$mom_bmi[which(pheno_clean$mom_bmi>100)]=NA


#mom's weight gain convert weight from lb to kg
pheno_clean$mom_wt_gain_kg=pheno_clean$wt_gain_select_2
pheno_clean$mom_wt_gain_kg[which(pheno_clean$mom_wt_gain_kg == -666 | pheno_clean$mom_wt_gain_kg == 888)]=NA
pheno_clean$wt_gain_select_1[which(pheno_clean$wt_gain_select_1 == -666 | pheno_clean$wt_gain_select_1 == 888)]=NA
tmp=pheno_clean$wt_gain_select_1* 0.45359237
pheno_clean$mom_wt_gain_kg[!is.na(tmp)]=tmp[!is.na(tmp)]
summary(pheno_clean$mom_wt_gain_kg)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.000   9.979  13.608  15.751  18.144  70.000     227 

#child's weight at birth, 1lb =16 ounces, 1lb=0.45359237kg=453.59237g,1oz=28.34952g
table(pheno_clean$child_bwt_select)
pheno_clean$child_bwt_g=pheno_clean$child_bwt_kg
pheno_clean$child_bwt_g[which(pheno_clean$child_bwt_g == -666 | pheno_clean$child_bwt_g == 888)]=NA
pheno_clean$child_bwt_select_1_lb[which(pheno_clean$child_bwt_select_1_lb == -666 | pheno_clean$child_bwt_select_1_lb == 888)]=NA
pheno_clean$child_bwt_select_1_oz[which(pheno_clean$child_bwt_select_1_oz == -666 | pheno_clean$child_bwt_select_1_oz == 888 | pheno_clean$child_bwt_select_1_oz == -999)]=NA
tmp=pheno_clean$child_bwt_select_1_lb * 453.59237 + pheno_clean$child_bwt_select_1_oz * 28.34952
pheno_clean$child_bwt_g[!is.na(tmp)]=tmp[!is.na(tmp)]
summary(pheno_clean$child_bwt_g)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 595.3  3005.0  3401.9  3367.6  3770.5  5358.1     127 

#to make things faster, let's select based on what we did before
pheno_clean_new=pheno_clean

pheno_clean=pheno_clean_new[,c(1,2,7,8,9,10,38:63,97:140)]
pheno_clean$alcohol_pregnancy[which(pheno_clean$alcohol_ever==0 & pheno_clean$alcohol_pregnancy==888)]=0
pheno_clean$alcohol_pregnancy[which(pheno_clean$alcohol_ever==0 & pheno_clean$alcohol_pregnancy== -666)]=0
pheno_clean$alcohol_frequency[which(pheno_clean$alcohol_ever==0 | pheno_clean$alcohol_pregnancy ==0)]=0


pheno_clean$smoke_preg[which(pheno_clean$smoke_ever==0 & pheno_clean$smoke_preg==888)]=0
pheno_clean$smoke_preg[which(pheno_clean$smoke_ever==0 & pheno_clean$smoke_preg== -666)]=0

#put eclampsia and preclam together
pheno_clean$eclampsia_preclam=pheno_clean$preclam
pheno_clean$eclampsia_preclam[which(pheno_clean$eclampsia==1)]=1

pheno_clean[pheno_clean == 999 | pheno_clean == 888 | pheno_clean == -999 | pheno_clean == -666] = NA

write.csv(pheno_clean,file="pheno_clean.csv")
save(pheno_clean,PGS_pheno,file="autism_model_data.RData")
