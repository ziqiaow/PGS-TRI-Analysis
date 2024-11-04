#ASD Data Analysis
#Data Cleaning for input for PGS-TRI
#April 24, 2024


setwd("./family/data/autism")
source("./family/R/PGScpt.R")

library(data.table)
library(fastDummies)

#load phenotype data
raw_data=read.csv("documents/alldata_July22.csv",stringsAsFactors = F)

#load calculated PRS based on 35,087 SNPs downloaded from PGS catalog
plink2=fread("plink/results/mergedscore.sscore")

#load family information
family=fread("plink/SPARK.iWES_v1.array.2022_02.fam",header = F)

#clean up family info
colnames(family)=c("FID","IID","father_id","mother_id","sex","phenotype")

#find complete trios
child_id=family$IID[which(family$father_id != "0" & family$mother_id != "0")]
any(duplicated(child_id))
any(duplicated(family$IID))
family_clean=family[match(child_id,family$IID),]
any(duplicated(family_clean$FID)) #there are related probands, i.e., from the same family
dim(family_clean)

#match with phenotype data
id=intersect(raw_data$Subject_ParticipantId,family_clean$IID)
family_clean=family_clean[match(id,family_clean$IID),]
dim(family_clean)

#create a PRS data frame
prs_pgscatalog=family_clean[,c(1,2,4,3,5)]
colnames(prs_pgscatalog)=c("FID","offspring_id","mother_id","father_id","sex")
prs_pgscatalog$offspring_prs=plink2$SCORE1_SUM[match(prs_pgscatalog$offspring_id,plink2$IID)]
prs_pgscatalog$mother_prs=plink2$SCORE1_SUM[match(prs_pgscatalog$mother_id,plink2$IID)]
prs_pgscatalog$father_prs=plink2$SCORE1_SUM[match(prs_pgscatalog$father_id,plink2$IID)]

prs_pgscatalog_complete=prs_pgscatalog[complete.cases(prs_pgscatalog),]
dim(prs_pgscatalog_complete)

#Clean up phenotype data
pheno=raw_data[match(prs_pgscatalog_complete$offspring_id,raw_data$Subject_ParticipantId),]
identical(pheno$Subject_ParticipantId,prs_pgscatalog_complete$offspring_id)
identical(pheno$Authorizer_ParticipantId,prs_pgscatalog_complete$mother_id)

#Add father's ID
pheno$father_id=prs_pgscatalog_complete$father_id

#Child's weight at birth, 1lb =16 ounces, 1lb=0.45359237kg=453.59237g, 1oz=28.34952g
table(pheno$child_bwt_select)
pheno$child_bwt_g=pheno$child_bwt_kg
pheno$child_bwt_g[which(pheno$child_bwt_g == -666 | pheno$child_bwt_g == 888)]=NA
pheno$child_bwt_select_1_lb[which(pheno$child_bwt_select_1_lb == -666 | pheno$child_bwt_select_1_lb == 888)]=NA
pheno$child_bwt_select_1_oz[which(pheno$child_bwt_select_1_oz == -666 | pheno$child_bwt_select_1_oz == 888 | pheno$child_bwt_select_1_oz == -999)]=NA
tmp=pheno$child_bwt_select_1_lb * 453.59237 + pheno$child_bwt_select_1_oz * 28.34952
pheno$child_bwt_g[!is.na(tmp)]=tmp[!is.na(tmp)]
summary(pheno$child_bwt_g)

pheno_clean=pheno[,c(1,5,3,97,7,8,9,10,39,41,42,43,44,45,46,48,49,50,52:63,98)]
pheno_clean$alcohol_pregnancy[which(pheno$alcohol_ever==0 & pheno$alcohol_pregnancy==888)]=0
pheno_clean$alcohol_pregnancy[which(pheno$alcohol_ever==0 & pheno$alcohol_pregnancy== -666)]=0
pheno_clean$alcohol_frequency[which(pheno_clean$alcohol_ever==0 | pheno_clean$alcohol_pregnancy ==0)]=0
pheno_clean$smoke_preg[which(pheno$smoke_ever==0 & pheno$smoke_preg==888)]=0
pheno_clean$smoke_preg[which(pheno$smoke_ever==0 & pheno$smoke_preg== -666)]=0

#Combine eclampsia and preclam
pheno_clean$eclampsia_preclam=pheno_clean$preclam
pheno_clean$eclampsia_preclam[which(pheno_clean$eclampsia==1)]=1

#Change codings to NA
pheno_clean[pheno_clean == 999 | pheno_clean == 888 | pheno_clean == -999 | pheno_clean == -666] = NA

#Combine phenotype data with genotype data for downstream analysis
identical(prs_pgscatalog_complete$offspring_id,pheno_clean$Subject_ParticipantId)
data.model=cbind(prs_pgscatalog_complete,pheno_clean)

#Remove children of the same family
dup=data.model[data.model$FID %in% unique(data.model$FID[ duplicated(data.model$FID)]) ,]
dup=dup[order(dup$FID),]
#Remove those with more NAs in phenotype data
id_rm1=c("SP0013622","SP0023047","SP0075872","SP0051559","SP0062146")
#Randomly generate 4 numbers of 1/2 to select to remove the rest related children
set.seed(04042023)
sample(1:2, 4, replace=TRUE)
id_rm2=c("SP0022322","SP0032771","SP0037276","SP0064631")
id_rm=c(id_rm1,id_rm2)
data.model=data.model[-match(id_rm,data.model$offspring_id),]


#combine smoke and tobacco
data.model$smoke_tobacco=data.model$smoke_ever
data.model$smoke_tobacco[which(data.model$tob_ever ==1 )]=1

data.model$smoke_tob_preg=data.model$smoke_preg
data.model$smoke_tob_preg[which(data.model$tob_preg ==1)]=1

#Combine depression and mental illness
data.model$depression_mentalill=data.model$depression
data.model$depression_mentalill[which(data.model$mental_illness==1)]=1

#Categorize the child weight to a binary variable given the small sample size (whether it is low birth weight or not)
data.model$child_bwt_low=data.model$child_bwt_g
data.model$child_bwt_low[which(data.model$child_bwt_g >= 2500)] = 0
data.model$child_bwt_low[which(data.model$child_bwt_g < 2500)] = 1

#create a new variable: alcohol: ever/preg: 0/0,1/0,1/1
#smoke_tobacco: ever/preg: 0/0,1/0,1/1
data.model$alcohol_ever_preg=data.model$alcohol_ever #group 0: 0/0
data.model$alcohol_ever_preg[which(data.model$alcohol_pregnancy == 0 & data.model$alcohol_ever == 1)] = 1
data.model$alcohol_ever_preg[which(data.model$alcohol_pregnancy == 1 & data.model$alcohol_ever == 1)] = 2

data.model$smoke_tob_ever_preg=data.model$smoke_tobacco #group 0: 0/0
data.model$smoke_tob_ever_preg[which(data.model$smoke_tob_preg == 0 & data.model$smoke_tobacco == 1)] = 1
data.model$smoke_tob_ever_preg[which(data.model$smoke_tob_preg == 1 & data.model$smoke_tobacco == 1)] = 2


#Load 1000 Genomes data
#load independent individual information for each race (children are removed)
eas=read.table("../1000G/1000G_EAS_ID.txt", header = F) 
eur=read.table("../1000G/1000G_EUR_ID.txt", header = F) 
sas=read.table("../1000G/1000G_SAS_ID.txt", header = F) 
afr=read.table("../1000G/1000G_AFR_ID.txt", header = F) 
amr=read.table("../1000G/1000G_AMR_ID.txt", header = F) 


plink2=fread("../1000G/autism/mergedscore.sscore")


#create a PRS data frame
prs_1000g_autism_EAS=plink2[match(eas$V1,plink2$IID),] #503
prs_1000g_autism_EUR=plink2[match(eur$V1,plink2$IID),] #498
prs_1000g_autism_SAS=plink2[match(sas$V1,plink2$IID),] #487
prs_1000g_autism_AFR=plink2[match(afr$V1,plink2$IID),] #659
prs_1000g_autism_AMR=plink2[match(amr$V1,plink2$IID),] #347

#Normalize by each race
data.model$offspring_prs_normalized_1000g_eas=( data.model$offspring_prs-mean(prs_1000g_autism_EAS$SCORE1_SUM) ) / sd(prs_1000g_autism_EAS$SCORE1_SUM)
data.model$mother_prs_normalized_1000g_eas=( data.model$mother_prs-mean(prs_1000g_autism_EAS$SCORE1_SUM) ) / sd(prs_1000g_autism_EAS$SCORE1_SUM)
data.model$father_prs_normalized_1000g_eas=( data.model$father_prs-mean(prs_1000g_autism_EAS$SCORE1_SUM) ) / sd(prs_1000g_autism_EAS$SCORE1_SUM)

data.model$offspring_prs_normalized_1000g_sas=( data.model$offspring_prs-mean(prs_1000g_autism_SAS$SCORE1_SUM) ) / sd(prs_1000g_autism_SAS$SCORE1_SUM)
data.model$mother_prs_normalized_1000g_sas=( data.model$mother_prs-mean(prs_1000g_autism_SAS$SCORE1_SUM) ) / sd(prs_1000g_autism_SAS$SCORE1_SUM)
data.model$father_prs_normalized_1000g_sas=( data.model$father_prs-mean(prs_1000g_autism_SAS$SCORE1_SUM) ) / sd(prs_1000g_autism_SAS$SCORE1_SUM)

data.model$offspring_prs_normalized_1000g_eur=( data.model$offspring_prs-mean(prs_1000g_autism_EUR$SCORE1_SUM) ) / sd(prs_1000g_autism_EUR$SCORE1_SUM)
data.model$mother_prs_normalized_1000g_eur=( data.model$mother_prs-mean(prs_1000g_autism_EUR$SCORE1_SUM) ) / sd(prs_1000g_autism_EUR$SCORE1_SUM)
data.model$father_prs_normalized_1000g_eur=( data.model$father_prs-mean(prs_1000g_autism_EUR$SCORE1_SUM) ) / sd(prs_1000g_autism_EUR$SCORE1_SUM)

data.model$offspring_prs_normalized_1000g_afr=( data.model$offspring_prs-mean(prs_1000g_autism_AFR$SCORE1_SUM) ) / sd(prs_1000g_autism_AFR$SCORE1_SUM)
data.model$mother_prs_normalized_1000g_afr=( data.model$mother_prs-mean(prs_1000g_autism_AFR$SCORE1_SUM) ) / sd(prs_1000g_autism_AFR$SCORE1_SUM)
data.model$father_prs_normalized_1000g_afr=( data.model$father_prs-mean(prs_1000g_autism_AFR$SCORE1_SUM) ) / sd(prs_1000g_autism_AFR$SCORE1_SUM)

data.model$offspring_prs_normalized_1000g_amr=( data.model$offspring_prs-mean(prs_1000g_autism_AMR$SCORE1_SUM) ) / sd(prs_1000g_autism_AMR$SCORE1_SUM)
data.model$mother_prs_normalized_1000g_amr=( data.model$mother_prs-mean(prs_1000g_autism_AMR$SCORE1_SUM) ) / sd(prs_1000g_autism_AMR$SCORE1_SUM)
data.model$father_prs_normalized_1000g_amr=( data.model$father_prs-mean(prs_1000g_autism_AMR$SCORE1_SUM) ) / sd(prs_1000g_autism_AMR$SCORE1_SUM)

all.prs.1000g=c(prs_1000g_autism_EAS$SCORE1_SUM,prs_1000g_autism_SAS$SCORE1_SUM,prs_1000g_autism_EUR$SCORE1_SUM,prs_1000g_autism_AFR$SCORE1_SUM,prs_1000g_autism_AMR$SCORE1_SUM)
data.model$offspring_prs_normalized_1000g=( data.model$offspring_prs-mean(all.prs.1000g) ) / sd(all.prs.1000g)
data.model$mother_prs_normalized_1000g=( data.model$mother_prs-mean(all.prs.1000g) ) / sd(all.prs.1000g)
data.model$father_prs_normalized_1000g=( data.model$father_prs-mean(all.prs.1000g) ) / sd(all.prs.1000g)

data.model$offspring_prs_normalized_1000g_combine=NA
data.model$mother_prs_normalized_1000g_combine=NA
data.model$father_prs_normalized_1000g_combine=NA
race=fread("documents/SPARK.iWES_v1.ancestry.2022_02.tsv",header = T)
race=race[match(data.model$offspring_id,race$spid),]
any(is.na(race$spid))
table(race$V3)
data.model$race=race$V3
data.model$race_detail=race$V2

#Remove the UNKNOWN ancestry individuals
data.model = data.model[-which(data.model$race=="UNKNOWN0" | data.model$race=="UNKNOWN1"),]

data.model$offspring_prs_normalized_1000g_combine[which(data.model$race=="EUR")]=data.model$offspring_prs_normalized_1000g_eur[which(data.model$race=="EUR")]
data.model$mother_prs_normalized_1000g_combine[which(data.model$race=="EUR")]=data.model$mother_prs_normalized_1000g_eur[which(data.model$race=="EUR")]
data.model$father_prs_normalized_1000g_combine[which(data.model$race=="EUR")]=data.model$father_prs_normalized_1000g_eur[which(data.model$race=="EUR")]
data.model$offspring_prs_normalized_1000g_combine[which(data.model$race=="EAS")]=data.model$offspring_prs_normalized_1000g_eas[which(data.model$race=="EAS")]
data.model$mother_prs_normalized_1000g_combine[which(data.model$race=="EAS")]=data.model$mother_prs_normalized_1000g_eas[which(data.model$race=="EAS")]
data.model$father_prs_normalized_1000g_combine[which(data.model$race=="EAS")]=data.model$father_prs_normalized_1000g_eas[which(data.model$race=="EAS")]
data.model$offspring_prs_normalized_1000g_combine[which(data.model$race=="SAS")]=data.model$offspring_prs_normalized_1000g_sas[which(data.model$race=="SAS")]
data.model$mother_prs_normalized_1000g_combine[which(data.model$race=="SAS")]=data.model$mother_prs_normalized_1000g_sas[which(data.model$race=="SAS")]
data.model$father_prs_normalized_1000g_combine[which(data.model$race=="SAS")]=data.model$father_prs_normalized_1000g_sas[which(data.model$race=="SAS")]
data.model$offspring_prs_normalized_1000g_combine[which(data.model$race=="AFR")]=data.model$offspring_prs_normalized_1000g_afr[which(data.model$race=="AFR")]
data.model$mother_prs_normalized_1000g_combine[which(data.model$race=="AFR")]=data.model$mother_prs_normalized_1000g_afr[which(data.model$race=="AFR")]
data.model$father_prs_normalized_1000g_combine[which(data.model$race=="AFR")]=data.model$father_prs_normalized_1000g_afr[which(data.model$race=="AFR")]
data.model$offspring_prs_normalized_1000g_combine[which(data.model$race=="AMR")]=data.model$offspring_prs_normalized_1000g_amr[which(data.model$race=="AMR")]
data.model$mother_prs_normalized_1000g_combine[which(data.model$race=="AMR")]=data.model$mother_prs_normalized_1000g_amr[which(data.model$race=="AMR")]
data.model$father_prs_normalized_1000g_combine[which(data.model$race=="AMR")]=data.model$father_prs_normalized_1000g_amr[which(data.model$race=="AMR")]


saveRDS(data.model,file="clean_data_asd_input.rds")
