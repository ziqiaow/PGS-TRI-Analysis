#ASD Data Analysis
#PRS-TRI for PGS main effect and parental nurturing effects
#April 24, 2024
#Ziqiao Wang


setwd("./family/data/autism")
source("/users/zwang4/family/R/PGS-TRI.R")

data.model <- readRDS("clean_data_asd_input.rds")


####################################################################
####################################################################
####################################################################
# Fit PRS using PRS-TRI, compared with pTDT test, for all ancestry groups combined
PRS.TRI.main=list()
PRS.TRI.main[[1]]=PGS.TRI(data.model$offspring_prs_normalized_1000g_combine,
                      data.model$mother_prs_normalized_1000g_combine,
                      data.model$father_prs_normalized_1000g_combine, 
                      GxE_int = FALSE, smalltriosize = T)
PRS.TRI.main[[1]]$res_beta
exp(0.2615448)

ptdt.test=list()
ptdt.test[[1]]=ptdt(data.model$offspring_prs_normalized_1000g_combine,
                    data.model$mother_prs_normalized_1000g_combine,
                    data.model$father_prs_normalized_1000g_combine)
ptdt.test[[1]]$res_beta
exp(0.1756665)



####################################################################
####################################################################
####################################################################
data.model1=data.model[,c(5,65:67,13,14,17,20,23,25,26,29,40,43:46,36,68)]


## European Ancestry
#First check PRS main effect
data.model2=data.model1[which(data.model1$race=="EUR"),]
PRS.TRI.main[[2]]=PGS.TRI(data.model2$offspring_prs_normalized_1000g_combine,
                      data.model2$mother_prs_normalized_1000g_combine,
                      data.model2$father_prs_normalized_1000g_combine, 
                      GxE_int = FALSE)
PRS.TRI.main[[2]]$res_beta
exp(0.2810934)

#check ptdt test
ptdt.test[[2]]=ptdt(data.model2$offspring_prs_normalized_1000g_combine,
                    data.model2$mother_prs_normalized_1000g_combine,
                    data.model2$father_prs_normalized_1000g_combine)
ptdt.test[[2]]$res_beta
exp(0.1895419)




################################################################
#################################################################
#################################################################

## African Ancestry
#First check PRS main effect
data.model2=data.model1[which(data.model1$race=="AFR"),]
PRS.TRI.main[[3]]=PGS.TRI(data.model2$offspring_prs_normalized_1000g_combine,
                      data.model2$mother_prs_normalized_1000g_combine,
                      data.model2$father_prs_normalized_1000g_combine, 
                      GxE_int = FALSE, smalltriosize = T)
PRS.TRI.main[[3]]$res_beta
exp(PRS.TRI.main[[3]]$res_beta[1])

#check ptdt test
ptdt.test[[3]]=ptdt(data.model2$offspring_prs_normalized_1000g_combine,
                    data.model2$mother_prs_normalized_1000g_combine,
                    data.model2$father_prs_normalized_1000g_combine)
ptdt.test[[3]]$res_beta
exp(ptdt.test[[3]]$res_beta[1])


################################################################
################################################################
################################################################
## Americas

#--------------------------------------
#First check PRS main effect
data.model2=data.model1[which(data.model1$race=="AMR"),]
PRS.TRI.main[[4]]=PGS.TRI(data.model2$offspring_prs_normalized_1000g_combine,
                      data.model2$mother_prs_normalized_1000g_combine,
                      data.model2$father_prs_normalized_1000g_combine, 
                      GxE_int = FALSE)
PRS.TRI.main[[4]]$res_beta
exp(PRS.TRI.main[[4]]$res_beta[1])

#check ptdt test
ptdt.test[[4]]=ptdt(data.model2$offspring_prs_normalized_1000g_combine,
                    data.model2$mother_prs_normalized_1000g_combine,
                    data.model2$father_prs_normalized_1000g_combine)
ptdt.test[[4]]$res_beta
exp(ptdt.test[[4]]$res_beta[1])



################################################################
################################################################
################################################################
## Asian Ancestry
#Combine South Asian and East Asian together due to small trio sizes

#--------------------------------------
#First check PRS main effect
data.model2=data.model1[which(data.model1$race=="SAS"| data.model1$race == "EAS"),]
PRS.TRI.main[[5]]=PGS.TRI(data.model2$offspring_prs_normalized_1000g_combine,
                      data.model2$mother_prs_normalized_1000g_combine,
                      data.model2$father_prs_normalized_1000g_combine, 
                      GxE_int = FALSE, smalltriosize = T)
PRS.TRI.main[[5]]$res_beta
exp(PRS.TRI.main[[5]]$res_beta[1])

#check ptdt test
ptdt.test[[5]]=ptdt(data.model2$offspring_prs_normalized_1000g_combine,
                    data.model2$mother_prs_normalized_1000g_combine,
                    data.model2$father_prs_normalized_1000g_combine)
ptdt.test[[5]]$res_beta
exp(ptdt.test[[5]]$res_beta[1])

# Remove AMR group in the analysis
data.model2=data.model1[which(data.model1$race != "AMR"),]
PRS.TRI.main[[6]]=PGS.TRI(data.model2$offspring_prs_normalized_1000g_combine,
                          data.model2$mother_prs_normalized_1000g_combine,
                          data.model2$father_prs_normalized_1000g_combine, 
                          GxE_int = FALSE)
PRS.TRI.main[[6]]$res_beta
exp(0.2918757)


ptdt.test[[6]]=ptdt(data.model2$offspring_prs_normalized_1000g_combine,
                    data.model2$mother_prs_normalized_1000g_combine,
                    data.model2$father_prs_normalized_1000g_combine)
ptdt.test[[6]]$res_beta
exp(0.1968166)



################################################################
################################################################
################################################################

# Check Parental Nurturing Effect
parent_nurture=list()
parent_nurture[[1]]=PGS.TRI(data.model$offspring_prs_normalized_1000g_combine,
                            data.model$mother_prs_normalized_1000g_combine,
                            data.model$father_prs_normalized_1000g_combine, 
                            parental_indirect = TRUE)
parent_nurture[[1]]

#Check different ancestry groups
# EUR
data.model1=data.model[which(data.model$race == "EUR"),]
parent_nurture[[2]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_combine,
                            data.model1$mother_prs_normalized_1000g_combine,
                            data.model1$father_prs_normalized_1000g_combine, 
                            parental_indirect = TRUE)
parent_nurture[[2]]

# AFR
data.model1=data.model[which(data.model$race == "AFR"),]
parent_nurture[[3]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_combine,
                            data.model1$mother_prs_normalized_1000g_combine,
                            data.model1$father_prs_normalized_1000g_combine, 
                            parental_indirect = TRUE, smalltriosize = T)
parent_nurture[[3]]

# AMR
data.model1=data.model[which(data.model$race == "AMR"),]
parent_nurture[[4]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_combine,
                            data.model1$mother_prs_normalized_1000g_combine,
                            data.model1$father_prs_normalized_1000g_combine, 
                            parental_indirect = TRUE)
parent_nurture[[4]]

# Asian (SAS+EAS)
data.model1=data.model[which(data.model$race == "SAS" | data.model$race == "EAS"),]
parent_nurture[[5]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_combine,
                            data.model1$mother_prs_normalized_1000g_combine,
                            data.model1$father_prs_normalized_1000g_combine, 
                            parental_indirect = TRUE, smalltriosize = T)
parent_nurture[[5]]

# All ancestry groups except AMR
data.model1=data.model[which(data.model$race != "AMR" ),]
parent_nurture[[6]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_combine,
                            data.model1$mother_prs_normalized_1000g_combine,
                            data.model1$father_prs_normalized_1000g_combine, 
                            parental_indirect = TRUE)
parent_nurture[[6]]












################################################
################################################
################################################
#Summarize the results
#PRS-TRI main PGS
tmp_res_main=array(0,c(6,3))
for(i in 1:6){
  tmp_res_main[i,1]=PRS.TRI.main[[i]]$res_beta[1]
  tmp_res_main[i,2]=PRS.TRI.main[[i]]$res_beta[2]
  tmp_res_main[i,3]=PRS.TRI.main[[i]]$res_beta[4]
  
}

rownames(tmp_res_main)=c("All Races", "EUR","AFR","AMR","SAS and EAS","All Races except AMR")
colnames(tmp_res_main)=c("beta","se","p")

write.csv(tmp_res_main,file="results/pgs_main_pgscpt.csv")
save(tmp_res_main,file="results/pgs_main_pgscpt.RData")

#ptdt main PGS
tmp_res_main=array(0,c(6,3))
for(i in 1:6){
  tmp_res_main[i,1]=ptdt.test[[i]]$res_beta[1]
  tmp_res_main[i,2]=ptdt.test[[i]]$res_beta[2]
  tmp_res_main[i,3]=ptdt.test[[i]]$res_beta[4]
  
}

rownames(tmp_res_main)=c("All Races", "EUR","AFR","AMR","SAS and EAS","All Races except AMR")
colnames(tmp_res_main)=c("beta","se","p")
write.csv(tmp_res_main,file="results/pgs_main_ptdt.csv")


#nurture effect
#summarize the main effect for our model
pgscpt_nurture=array(0,c(6,3))
for(i in 1:6){
  pgscpt_nurture[i,]=parent_nurture[[i]]$res_delta[c(1,2,4)]
  
  
}
rownames(pgscpt_nurture)=c("All Races", "EUR","AFR","AMR","SAS and EAS","All Races except AMR")
colnames(pgscpt_nurture)=c("beta","se","p")

save(pgscpt_nurture,file="results/nurture.RData")
write.csv(pgscpt_nurture,file="results/nurture.csv")
