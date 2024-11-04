#---------------------------------------------
#summarize the results for OC
#August 16, 2024
#---------------------------------------------

setwd("./family/data/oral_cleft")
load("data_clean_1000g.RData")
source("/users/zwang4/family/R/PGS-TRI.R")

#remove unknown subtype of OC
data.model=data.model[-which(data.model$proband_clp=="9"),]
# Check OR for PRS alone using our developed likelihood-based method, compared with pTDT test
ptdt.test.all=list()
pgscpt.test.all=list()
ptdt.test.eur=list()
pgscpt.test.eur=list()
ptdt.test.eas=list()
pgscpt.test.eas=list()
#All (EUR and ASA)
pgscpt.test.all[[1]]=PGS.TRI(data.model$offspring_prs_normalized_1000g_eas_eur_combine,
                             data.model$mother_prs_normalized_1000g_eas_eur_combine,
                             data.model$father_prs_normalized_1000g_eas_eur_combine, 
                             GxE_int = FALSE)
pgscpt.test.all[[1]]$res_beta
exp(0.5012097)

#Check the method in Weiner et al., Nat Genet 2017
ptdt.test.all[[1]]=ptdt(data.model$offspring_prs_normalized_1000g_eas_eur_combine,
                        data.model$mother_prs_normalized_1000g_eas_eur_combine,
                        data.model$father_prs_normalized_1000g_eas_eur_combine)
ptdt.test.all[[1]]$res_beta
exp(0.3590803)

#All, EUR
pgscpt.test.eur[[1]]=PGS.TRI(data.model$offspring_prs_normalized_1000g_eas_eur_combine[which(data.model$proband_race2=="EUR")],
                             data.model$mother_prs_normalized_1000g_eas_eur_combine[which(data.model$proband_race2=="EUR")],
                             data.model$father_prs_normalized_1000g_eas_eur_combine[which(data.model$proband_race2=="EUR")], 
                             GxE_int = FALSE)

pgscpt.test.eur[[1]]$res_beta
#Check the method in Weiner et al., Nat Genet 2017
(ptdt.test.eur[[1]]=ptdt(data.model$offspring_prs_normalized_1000g_eas_eur_combine[which(data.model$proband_race2=="EUR")],
                         data.model$mother_prs_normalized_1000g_eas_eur_combine[which(data.model$proband_race2=="EUR")],
                         data.model$father_prs_normalized_1000g_eas_eur_combine[which(data.model$proband_race2=="EUR")]))



#All, EAS
pgscpt.test.eas[[1]]=PGS.TRI(data.model$offspring_prs_normalized_1000g_eas_eur_combine[which(data.model$proband_race2=="ASA")],
                             data.model$mother_prs_normalized_1000g_eas_eur_combine[which(data.model$proband_race2=="ASA")],
                             data.model$father_prs_normalized_1000g_eas_eur_combine[which(data.model$proband_race2=="ASA")], 
                             GxE_int = FALSE)
pgscpt.test.eas[[1]]$res_beta

#Check the method in Weiner et al., Nat Genet 2017
(ptdt.test.eas[[1]]=ptdt(data.model$offspring_prs_normalized_1000g_eas_eur_combine[which(data.model$proband_race2=="ASA")],
                         data.model$mother_prs_normalized_1000g_eas_eur_combine[which(data.model$proband_race2=="ASA")],
                         data.model$father_prs_normalized_1000g_eas_eur_combine[which(data.model$proband_race2=="ASA")]))



#Check the subtypes
#1=Cleft Palate only	2=Cleft Lip only	3=Cleft Lip & Palate	9=Unknown
table(data.model$proband_clp)

#Group cleft lip only and cleft lip & palate together to CL/P 
data.model$clp=data.model$proband_clp
data.model$clp[which(data.model$proband_clp == "2" | data.model$proband_clp == "3")]="CL/P"
data.model$clp[which(data.model$proband_clp=="1")]="CP"
data.model$clp[which(data.model$proband_clp=="9")]=NA
table(data.model$clp) #CP: cleft palate only; CL/P: cleft lip with/without palate
table(data.model$clp,data.model$proband_race2)
#Check for each race and disease subtype
#EUR/CP
data.model1=data.model[which(data.model$proband_race2=="EUR" & data.model$clp=="CP"),]
pgscpt.test.eur[[5]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                             data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                             data.model1$father_prs_normalized_1000g_eas_eur_combine, 
                             GxE_int = FALSE)
pgscpt.test.eur[[5]]$res_beta
#Check the method in Weiner et al., Nat Genet 2017
(ptdt.test.eur[[5]]=ptdt(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                         data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                         data.model1$father_prs_normalized_1000g_eas_eur_combine))



#EUR/CL/P
data.model1=data.model[which(data.model$proband_race2=="EUR" & data.model$clp=="CL/P"),]
pgscpt.test.eur[[2]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                             data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                             data.model1$father_prs_normalized_1000g_eas_eur_combine, 
                             GxE_int = FALSE)
pgscpt.test.eur[[2]]$res_beta
(ptdt.test.eur[[2]]=ptdt(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                         data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                         data.model1$father_prs_normalized_1000g_eas_eur_combine))


#EUR/CL only
data.model1=data.model[which(data.model$proband_race2=="EUR" & data.model$proband_clp=="2"),]
pgscpt.test.eur[[3]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                             data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                             data.model1$father_prs_normalized_1000g_eas_eur_combine, 
                             GxE_int = FALSE)
pgscpt.test.eur[[3]]$res_beta
(ptdt.test.eur[[3]]=ptdt(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                         data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                         data.model1$father_prs_normalized_1000g_eas_eur_combine))


#EUR/CL&P 
data.model1=data.model[which(data.model$proband_race2=="EUR" & data.model$proband_clp=="3"),]
pgscpt.test.eur[[4]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                             data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                             data.model1$father_prs_normalized_1000g_eas_eur_combine, 
                             GxE_int = FALSE)
pgscpt.test.eur[[4]]$res_beta
(ptdt.test.eur[[4]]=ptdt(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                         data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                         data.model1$father_prs_normalized_1000g_eas_eur_combine))


#ASA/CP
data.model1=data.model[which(data.model$proband_race2=="ASA" & data.model$clp=="CP"),]
pgscpt.test.eas[[5]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                             data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                             data.model1$father_prs_normalized_1000g_eas_eur_combine, 
                             GxE_int = FALSE)
pgscpt.test.eas[[5]]$res_beta
(ptdt.test.eas[[5]]=ptdt(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                         data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                         data.model1$father_prs_normalized_1000g_eas_eur_combine))


#ASA/CL/P
data.model1=data.model[which(data.model$proband_race2=="ASA" & data.model$clp=="CL/P"),]
pgscpt.test.eas[[2]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                             data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                             data.model1$father_prs_normalized_1000g_eas_eur_combine, 
                             GxE_int = FALSE)
pgscpt.test.eas[[2]]$res_beta
(ptdt.test.eas[[2]]=ptdt(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                         data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                         data.model1$father_prs_normalized_1000g_eas_eur_combine))


#ASA/CL only
data.model1=data.model[which(data.model$proband_race2=="ASA" & data.model$proband_clp=="2"),]
pgscpt.test.eas[[3]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                             data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                             data.model1$father_prs_normalized_1000g_eas_eur_combine, 
                             GxE_int = FALSE)
pgscpt.test.eas[[3]]$res_beta
(ptdt.test.eas[[3]]=ptdt(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                         data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                         data.model1$father_prs_normalized_1000g_eas_eur_combine))


#ASA/CL&P 
data.model1=data.model[which(data.model$proband_race2=="ASA" & data.model$proband_clp=="3"),]
pgscpt.test.eas[[4]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                             data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                             data.model1$father_prs_normalized_1000g_eas_eur_combine, 
                             GxE_int = FALSE)
pgscpt.test.eas[[4]]$res_beta
(ptdt.test.eas[[4]]=ptdt(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                         data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                         data.model1$father_prs_normalized_1000g_eas_eur_combine))


#ALL/CP
data.model1=data.model[which( data.model$clp=="CP"),]
pgscpt.test.all[[5]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                             data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                             data.model1$father_prs_normalized_1000g_eas_eur_combine, 
                             GxE_int = FALSE)
pgscpt.test.all[[5]]$res_beta
(ptdt.test.all[[5]]=ptdt(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                         data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                         data.model1$father_prs_normalized_1000g_eas_eur_combine))


#ALL/CL/P
data.model1=data.model[which(data.model$clp=="CL/P"),]
pgscpt.test.all[[2]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                             data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                             data.model1$father_prs_normalized_1000g_eas_eur_combine, 
                             GxE_int = FALSE)
pgscpt.test.all[[2]]$res_beta
(ptdt.test.all[[2]]=ptdt(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                         data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                         data.model1$father_prs_normalized_1000g_eas_eur_combine))


#ALL/CL only
data.model1=data.model[which(data.model$proband_clp=="2"),]
pgscpt.test.all[[3]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                             data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                             data.model1$father_prs_normalized_1000g_eas_eur_combine, 
                             GxE_int = FALSE)
pgscpt.test.all[[3]]$res_beta
(ptdt.test.all[[3]]=ptdt(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                         data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                         data.model1$father_prs_normalized_1000g_eas_eur_combine))


#ALL/CL&P 
data.model1=data.model[which(data.model$proband_clp=="3"),]
pgscpt.test.all[[4]]=PGS.TRI(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                             data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                             data.model1$father_prs_normalized_1000g_eas_eur_combine, 
                             GxE_int = FALSE)
pgscpt.test.all[[4]]$res_beta
(ptdt.test.all[[4]]=ptdt(data.model1$offspring_prs_normalized_1000g_eas_eur_combine,
                         data.model1$mother_prs_normalized_1000g_eas_eur_combine,
                         data.model1$father_prs_normalized_1000g_eas_eur_combine))

#summarize the main effect for our model
pgscpt_main=array(0,c(5,9))
for(i in 1:5){
  pgscpt_main[i,1:3]=pgscpt.test.all[[i]]$res_beta[c(1,2,4)]
  pgscpt_main[i,4:6]=pgscpt.test.eur[[i]]$res_beta[c(1,2,4)]
  pgscpt_main[i,7:9]=pgscpt.test.eas[[i]]$res_beta[c(1,2,4)]
  
}
rownames(pgscpt_main)=c("all","CL/P","CL","CL&P","CP")
colnames(pgscpt_main)=c("all_beta","all_se","all_p",
                        "eur_beta","eur_se","eur_p",
                        "asian_beta","asian_se","asian_p")


ptdt_main=array(0,c(5,9))
for(i in 1:5){
  ptdt_main[i,1:3]=ptdt.test.all[[i]]$res_beta[c(1,2,4)]
  ptdt_main[i,4:6]=ptdt.test.eur[[i]]$res_beta[c(1,2,4)]
  ptdt_main[i,7:9]=ptdt.test.eas[[i]]$res_beta[c(1,2,4)]
  
}
rownames(ptdt_main)=c("all","CL/P","CL","CL&P","CP")
colnames(ptdt_main)=c("all_beta","all_se","all_p",
                        "eur_beta","eur_se","eur_p",
                        "asian_beta","asian_se","asian_p")

save(ptdt_main,pgscpt_main,file="results/PGS_main_cpt_ptdt.RData")

# Check Parental Nurturing Effect
pgscpt.test.all=list()
pgscpt.test.eur=list()
pgscpt.test.eas=list()

pgscpt.test.all[[1]]=PGS.TRI(pgs_offspring = data.model$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)

#EUR
data.model1=data.model[which(data.model$proband_race2=="EUR"),]
pgscpt.test.eur[[1]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)

#EAS
data.model1=data.model[which(data.model$proband_race2=="ASA"),]
pgscpt.test.eas[[1]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)


#Check different subtype and ancestry group combinations
#CP, ALL
data.model1=data.model[which(data.model$proband_clp=="1"),]
pgscpt.test.all[[5]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)

#CL/P (with or without palate), ALL
data.model1=data.model[which(data.model$clp=="CL/P"),]
pgscpt.test.all[[2]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)

#CL alone, ALL
data.model1=data.model[which(data.model$proband_clp=="2"),]
pgscpt.test.all[[3]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)

#CL and palate, ALL
data.model1=data.model[which(data.model$proband_clp=="3"),]
pgscpt.test.all[[4]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)

#CP, EUR
data.model1=data.model[which(data.model$proband_clp=="1" & data.model$proband_race2=="EUR"),]
pgscpt.test.eur[[5]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)

#CL/P (with or without palate), EUR
data.model1=data.model[which(data.model$clp=="CL/P" & data.model$proband_race2=="EUR"),]
pgscpt.test.eur[[2]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)

#CL alone, EUR
data.model1=data.model[which(data.model$proband_clp=="2" & data.model$proband_race2=="EUR"),]
pgscpt.test.eur[[3]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)

#CL and palate, EUR
data.model1=data.model[which(data.model$proband_clp=="3" & data.model$proband_race2=="EUR"),]
pgscpt.test.eur[[4]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)
#CP, EAS
data.model1=data.model[which(data.model$proband_clp=="1" & data.model$proband_race2=="ASA"),]
pgscpt.test.eas[[5]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)

#CL/P (with or without palate), EAS
data.model1=data.model[which(data.model$clp=="CL/P" & data.model$proband_race2=="ASA"),]
pgscpt.test.eas[[2]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)

#CL alone, EAS
data.model1=data.model[which(data.model$proband_clp=="2" & data.model$proband_race2=="ASA"),]
pgscpt.test.eas[[3]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)


#CL and palate, EAS
data.model1=data.model[which(data.model$proband_clp=="3" & data.model$proband_race2=="ASA"),]
pgscpt.test.eas[[4]]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
        pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
        pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
        parental_indirect = TRUE)


#summarize the main effect for our model
pgscpt_nurture=array(0,c(5,9))
for(i in 1:5){
  pgscpt_nurture[i,1:3]=pgscpt.test.all[[i]]$res_delta[c(1,2,4)]
  pgscpt_nurture[i,4:6]=pgscpt.test.eur[[i]]$res_delta[c(1,2,4)]
  pgscpt_nurture[i,7:9]=pgscpt.test.eas[[i]]$res_delta[c(1,2,4)]
  
}
rownames(pgscpt_nurture)=c("all","CL/P","CL","CL&P","CP")
colnames(pgscpt_nurture)=c("all_beta","all_se","all_p",
                        "eur_beta","eur_se","eur_p",
                        "asian_beta","asian_se","asian_p")

save(pgscpt_nurture,file="results/nurture_ptdt.RData")
write.csv(pgscpt_nurture,file="results/nurture_ptdt.csv")


# PRS x Maternal Exposures using our developed model

res_all_CLP=list()
for(i in 1:3){
  tmp=array(0,c(4,3))
  data.model1=data.model[which(data.model$proband_clp==i),]
  data.model1[data.model1 == "9"]=NA
  
  tmp[1,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
          pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
          pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
          GxE_int = TRUE,
          formula = ~ msmoke,
          E = data.model1, 
          side = 2, 
          numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp[2,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
          pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
          pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
          GxE_int = TRUE,
          formula = ~ ets,
          E = data.model1, 
          side = 2, 
          numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp[3,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
          pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
          pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
          GxE_int = TRUE,
          formula = ~ mvits,
          E = data.model1, 
          side = 2, 
          numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp[4,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
          pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
          pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
          GxE_int = TRUE,
          formula = ~ mdrink,
          E = data.model1, 
          side = 2, 
          numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eur=array(0,c(4,3))
  data.model1=data.model[which(data.model$proband_clp==i  & data.model$proband_race2=="EUR"),]
  data.model1[data.model1 == "9"]=NA
  
  tmp_eur[1,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                  pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                  pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                  GxE_int = TRUE,
                  formula = ~ msmoke,
                  E = data.model1, 
                  side = 2, 
                  numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eur[2,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                  pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                  pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                  GxE_int = TRUE,
                  formula = ~ ets,
                  E = data.model1, 
                  side = 2, 
                  numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eur[3,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                  pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                  pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                  GxE_int = TRUE,
                  formula = ~ mvits,
                  E = data.model1, 
                  side = 2, 
                  numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eur[4,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                  pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                  pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                  GxE_int = TRUE,
                  formula = ~ mdrink,
                  E = data.model1, 
                  side = 2, 
                  numDeriv = F)$res_beta[2,c(1,2,4)]
  
  
  tmp_eas=array(0,c(4,3))
  data.model1=data.model[which(data.model$proband_clp==i  & data.model$proband_race2=="ASA"),]
  data.model1[data.model1 == "9"]=NA
  
  tmp_eas[1,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ msmoke,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eas[2,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ ets,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eas[3,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ mvits,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eas[4,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ mdrink,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_all=cbind(tmp,tmp_eur,tmp_eas)
  colnames(tmp_all)=c("all_beta","all_se","all_p","eur_beta","eur_se","eur_p","eas_beta","eas_se","eas_p")
  rownames(tmp_all)=c("smoke","ets","vit","drink")
  res_all_CLP[[i]]=tmp_all
}

res_all_CLP[[2]][which(res_all_CLP[[2]]==0)]=NA


i=4
  tmp=array(0,c(4,3))
  data.model1=data.model[which(data.model$clp=="CL/P"),]
  data.model1[data.model1 == "9"]=NA
  
  tmp[1,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                  pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                  pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                  GxE_int = TRUE,
                  formula = ~ msmoke,
                  E = data.model1, 
                  side = 2, 
                  numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp[2,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                  pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                  pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                  GxE_int = TRUE,
                  formula = ~ ets,
                  E = data.model1, 
                  side = 2, 
                  numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp[3,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                  pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                  pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                  GxE_int = TRUE,
                  formula = ~ mvits,
                  E = data.model1, 
                  side = 2, 
                  numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp[4,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                  pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                  pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                  GxE_int = TRUE,
                  formula = ~ mdrink,
                  E = data.model1, 
                  side = 2, 
                  numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eur=array(0,c(4,3))
  data.model1=data.model[which(data.model$clp=="CL/P"  & data.model$proband_race2=="EUR"),]
  data.model1[data.model1 == "9"]=NA
  
  tmp_eur[1,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ msmoke,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eur[2,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ ets,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eur[3,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ mvits,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eur[4,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ mdrink,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  
  tmp_eas=array(0,c(4,3))
  data.model1=data.model[which(data.model$clp=="CL/P"  & data.model$proband_race2=="ASA"),]
  data.model1[data.model1 == "9"]=NA
  
  tmp_eas[1,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ msmoke,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eas[2,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ ets,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eas[3,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ mvits,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eas[4,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ mdrink,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_all=cbind(tmp,tmp_eur,tmp_eas)
  colnames(tmp_all)=c("all_beta","all_se","all_p","eur_beta","eur_se","eur_p","eas_beta","eas_se","eas_p")
rownames(tmp_all)=c("smoke","ets","vit","drink")
  res_all_CLP[[i]]=tmp_all

  
  i=5
  tmp=array(0,c(4,3))
  data.model1=data.model
  data.model1[data.model1 == "9"]=NA
  
  tmp[1,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                  pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                  pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                  GxE_int = TRUE,
                  formula = ~ msmoke,
                  E = data.model1, 
                  side = 2, 
                  numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp[2,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                  pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                  pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                  GxE_int = TRUE,
                  formula = ~ ets,
                  E = data.model1, 
                  side = 2, 
                  numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp[3,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                  pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                  pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                  GxE_int = TRUE,
                  formula = ~ mvits,
                  E = data.model1, 
                  side = 2, 
                  numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp[4,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                  pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                  pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                  GxE_int = TRUE,
                  formula = ~ mdrink,
                  E = data.model1, 
                  side = 2, 
                  numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eur=array(0,c(4,3))
  data.model1=data.model[which(data.model$proband_race2=="EUR"),]
  data.model1[data.model1 == "9"]=NA
  
  tmp_eur[1,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ msmoke,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eur[2,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ ets,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eur[3,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ mvits,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eur[4,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ mdrink,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  
  tmp_eas=array(0,c(4,3))
  data.model1=data.model[which(data.model$proband_race2=="ASA"),]
  data.model1[data.model1 == "9"]=NA
  
  tmp_eas[1,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ msmoke,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eas[2,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ ets,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eas[3,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ mvits,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_eas[4,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                      pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                      pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                      GxE_int = TRUE,
                      formula = ~ mdrink,
                      E = data.model1, 
                      side = 2, 
                      numDeriv = F)$res_beta[2,c(1,2,4)]
  
  tmp_all=cbind(tmp,tmp_eur,tmp_eas)
  colnames(tmp_all)=c("all_beta","all_se","all_p","eur_beta","eur_se","eur_p","eas_beta","eas_se","eas_p")
  rownames(tmp_all)=c("smoke","ets","vit","drink")
  res_all_CLP[[i]]=tmp_all
  
  
names(res_all_CLP)=c("CP","CL","CL&P","CL/P","all")
res_interaction=data.frame(do.call("rbind", res_all_CLP))
res_interaction$subtype=c(rep(names(res_all_CLP),each=4))
res_interaction$envir=c(rep(rownames(res_interaction)[1:4],5))
save(res_interaction,file="results/interaction.RData")

res_interaction$envir=factor(res_interaction$envir,levels=c("smoke","ets","vit","drink"))
res_interaction$subtype=factor(res_interaction$subtype,levels=c("CL/P","CL","CL&P","all","CP"))
res_interaction=res_interaction[order(res_interaction$envir,res_interaction$subtype),]
write.csv(res_interaction,file="results/interaction.csv")

#For CL/P, check sex interactions
tmp=array(0,c(5,3))
data.model1=data.model[which(data.model$clp=="CL/P"),]
data.model1[data.model1 == "9"]=NA

tmp[1,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                GxE_int = TRUE,
                formula = ~ msmoke,
                E = data.model1, 
                side = 2, 
                numDeriv = F)$res_beta[2,c(1,2,4)]

tmp[2,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                GxE_int = TRUE,
                formula = ~ ets,
                E = data.model1, 
                side = 2, 
                numDeriv = F)$res_beta[2,c(1,2,4)]

tmp[3,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                GxE_int = TRUE,
                formula = ~ mvits,
                E = data.model1, 
                side = 2, 
                numDeriv = F)$res_beta[2,c(1,2,4)]

tmp[4,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                GxE_int = TRUE,
                formula = ~ mdrink,
                E = data.model1, 
                side = 2, 
                numDeriv = F)$res_beta[2,c(1,2,4)]

tmp[5,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                GxE_int = TRUE,
                formula = ~ factor(sex),
                E = data.model1, 
                side = 2, 
                numDeriv = F)$res_beta[2,c(1,2,4)]

tmp_eur=array(0,c(5,3))
data.model1=data.model[which(data.model$clp=="CL/P"  & data.model$proband_race2=="EUR"),]
data.model1[data.model1 == "9"]=NA

tmp_eur[1,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                    pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                    pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                    GxE_int = TRUE,
                    formula = ~ msmoke,
                    E = data.model1, 
                    side = 2, 
                    numDeriv = F)$res_beta[2,c(1,2,4)]

tmp_eur[2,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                    pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                    pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                    GxE_int = TRUE,
                    formula = ~ ets,
                    E = data.model1, 
                    side = 2, 
                    numDeriv = F)$res_beta[2,c(1,2,4)]

tmp_eur[3,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                    pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                    pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                    GxE_int = TRUE,
                    formula = ~ mvits,
                    E = data.model1, 
                    side = 2, 
                    numDeriv = F)$res_beta[2,c(1,2,4)]

tmp_eur[4,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                    pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                    pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                    GxE_int = TRUE,
                    formula = ~ mdrink,
                    E = data.model1, 
                    side = 2, 
                    numDeriv = F)$res_beta[2,c(1,2,4)]


tmp_eur[5,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                    pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                    pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                    GxE_int = TRUE,
                    formula = ~ factor(sex),
                    E = data.model1, 
                    side = 2, 
                    numDeriv = F)$res_beta[2,c(1,2,4)]


tmp_eas=array(0,c(5,3))
data.model1=data.model[which(data.model$clp=="CL/P"  & data.model$proband_race2=="ASA"),]
data.model1[data.model1 == "9"]=NA

tmp_eas[1,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                    pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                    pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                    GxE_int = TRUE,
                    formula = ~ msmoke,
                    E = data.model1, 
                    side = 2, 
                    numDeriv = F)$res_beta[2,c(1,2,4)]

tmp_eas[2,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                    pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                    pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                    GxE_int = TRUE,
                    formula = ~ ets,
                    E = data.model1, 
                    side = 2, 
                    numDeriv = F)$res_beta[2,c(1,2,4)]

tmp_eas[3,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                    pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                    pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                    GxE_int = TRUE,
                    formula = ~ mvits,
                    E = data.model1, 
                    side = 2, 
                    numDeriv = F)$res_beta[2,c(1,2,4)]

tmp_eas[4,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                    pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                    pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                    GxE_int = TRUE,
                    formula = ~ mdrink,
                    E = data.model1, 
                    side = 2, 
                    numDeriv = F)$res_beta[2,c(1,2,4)]

tmp_eas[5,]=PGS.TRI(pgs_offspring = data.model1$offspring_prs_normalized_1000g_eas_eur_combine, 
                    pgs_mother = data.model1$mother_prs_normalized_1000g_eas_eur_combine, 
                    pgs_father = data.model1$father_prs_normalized_1000g_eas_eur_combine,
                    GxE_int = TRUE,
                    formula = ~ factor(sex),
                    E = data.model1, 
                    side = 2, 
                    numDeriv = F)$res_beta[2,c(1,2,4)]

tmp_all=cbind(tmp,tmp_eur,tmp_eas)
colnames(tmp_all)=c("all_beta","all_se","all_p","eur_beta","eur_se","eur_p","eas_beta","eas_se","eas_p")
rownames(tmp_all)=c("smoke","ets","vit","drink","sex")
interaction_CLP=tmp_all
save(interaction_CLP,file="results/interaction_clp.RData")
write.csv(interaction_CLP,file="results/interaction_clp.csv")
