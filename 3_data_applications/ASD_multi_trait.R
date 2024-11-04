#summary statistics can be found at PGS catalog
#or
#https://pgc.unc.edu/for-researchers/download-results/
setwd("./family/data/pleiotropy")
trait=c("edu","schizophrenia","depression","bipolar","adhd","bmi","bmi_prive","bipolar1","bipolar2","neuroticism","insomnia","chronotype")

results_all=array(0,c(length(trait),19))
results_all=data.frame(results_all)
rownames(results_all)=trait
colnames(results_all)=c("All.noAMR.beta","All.noAMR.se","All.noAMR.pvalue",
                        "EUR.beta","EUR.se","EUR.pvalue",
                        "AFR.beta","AFR.se","AFR.pvalue",
                        "AMR.beta","AMR.se","AMR.pvalue",
                        "AS.beta","AS.se","AS.pvalue",
                        "All.beta","All.se","All.pvalue",
                        "#SNPs")
source("./family/R/PRS-TRI.R")
#N (All) = 1364; N (EUR) = 1250

delta_all=array(0,c(length(trait),19))
delta_all=data.frame(delta_all)
rownames(delta_all)=trait

colnames(delta_all)=c("All.noAMR.beta","All.noAMR.se","All.noAMR.pvalue",
                      "EUR.beta","EUR.se","EUR.pvalue",
                      "AFR.beta","AFR.se","AFR.pvalue",
                      "AMR.beta","AMR.se","AMR.pvalue",
                      "AS.beta","AS.se","AS.pvalue",
                      "All.beta","All.se","All.pvalue",
                      "#SNPs")


h=0
for(i in trait){
  if(i != "adhd"){
    load(paste0("./family/data/pleiotropy/",i,"/",i,"_clean_trio_standardized.RData"))
  } else {load(paste0("./family/data/pleiotropy/adhd_redo/PRScs/adhd_trio_standardized_prscs.RData"))
  }  
  h=h+1
       data.model2=data.model1[-which(data.model1$race=="AMR"),]
       pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                        data.model2$mother_prs_standardized,
                        data.model2$father_prs_standardized, 
                        GxE_int = FALSE,
                        parental_nurture = T)
       #results_all[h,1:3]=pgs_main$res_beta[c(1,2,4)]
       delta_all[h,1:3]=pgs_main$res_delta[c(1,2,4)]
       
       #EUR
       data.model2=data.model1[which(data.model1$race=="EUR"),]
       pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                        data.model2$mother_prs_standardized,
                        data.model2$father_prs_standardized, 
                        GxE_int = FALSE,
                        parental_nurture = T)
       #results_all[h,4:6]=pgs_main$res_beta[c(1,2,4)]
       delta_all[h,4:6]=pgs_main$res_delta[c(1,2,4)]
       
       
       #AFR
       data.model2=data.model1[which(data.model1$race=="AFR"),]
       pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                        data.model2$mother_prs_standardized,
                        data.model2$father_prs_standardized, 
                        GxE_int = FALSE,
                        parental_nurture = T)
       #results_all[h,7:9]=pgs_main$res_beta[c(1,2,4)]
       delta_all[h,7:9]=pgs_main$res_delta[c(1,2,4)]
       
       #AMR
       data.model2=data.model1[which(data.model1$race=="AMR"),]
       pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                        data.model2$mother_prs_standardized,
                        data.model2$father_prs_standardized, 
                        GxE_int = FALSE,
                        parental_nurture = T)
       #results_all[h,10:12]=pgs_main$res_beta[c(1,2,4)]
       delta_all[h,10:12]=pgs_main$res_delta[c(1,2,4)]
       
       #AS
       data.model2=data.model1[which(data.model1$race=="SAS" | data.model1$race=="EAS"),]
       pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                        data.model2$mother_prs_standardized,
                        data.model2$father_prs_standardized, 
                        GxE_int = FALSE,
                        parental_nurture = T)
       #results_all[h,13:15]=pgs_main$res_beta[c(1,2,4)]
       delta_all[h,13:15]=pgs_main$res_delta[c(1,2,4)]
       
       
       #all
       data.model2=data.model1
       pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                        data.model2$mother_prs_standardized,
                        data.model2$father_prs_standardized, 
                        GxE_int = FALSE,
                        parental_nurture = T)
       #results_all[h,16:18]=pgs_main$res_beta[c(1,2,4)]
       delta_all[h,16:18]=pgs_main$res_delta[c(1,2,4)]
       
}



h=0
for(i in trait){
  if(i != "adhd"){
    load(paste0("./family/data/pleiotropy/",i,"/",i,"_clean_trio_standardized.RData"))
  } else {load(paste0("./family/data/pleiotropy/adhd_redo/PRScs/adhd_trio_standardized_prscs.RData"))
  }  
  h=h+1
  data.model2=data.model1[-which(data.model1$race=="AMR"),]
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_nurture = F)
  results_all[h,1:3]=pgs_main$res_beta[c(1,2,4)]
 # delta_all[h,1:3]=pgs_main$res_delta[c(1,2,4)]
  
  #EUR
  data.model2=data.model1[which(data.model1$race=="EUR"),]
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_nurture = F)
  results_all[h,4:6]=pgs_main$res_beta[c(1,2,4)]
 # delta_all[h,4:6]=pgs_main$res_delta[c(1,2,4)]
  
  
  #AFR
  data.model2=data.model1[which(data.model1$race=="AFR"),]
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_nurture = F)
  results_all[h,7:9]=pgs_main$res_beta[c(1,2,4)]
 # delta_all[h,7:9]=pgs_main$res_delta[c(1,2,4)]
  
  #AMR
  data.model2=data.model1[which(data.model1$race=="AMR"),]
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_nurture = F)
  results_all[h,10:12]=pgs_main$res_beta[c(1,2,4)]
 # delta_all[h,10:12]=pgs_main$res_delta[c(1,2,4)]
  
  #AS
  data.model2=data.model1[which(data.model1$race=="SAS" | data.model1$race=="EAS"),]
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_nurture = F)
  results_all[h,13:15]=pgs_main$res_beta[c(1,2,4)]
  #delta_all[h,13:15]=pgs_main$res_delta[c(1,2,4)]
  
  
  #all
  data.model2=data.model1
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_nurture = F)
  results_all[h,16:18]=pgs_main$res_beta[c(1,2,4)]
 # delta_all[h,16:18]=pgs_main$res_delta[c(1,2,4)]
  
}
results_all$`#SNPs`=delta_all$`#SNPs`=c(49265,585938,18746,741261,419357,1414,107018,  731828, 730189 ,53398,900119,927726)
results_all
delta_all
save(results_all,delta_all,file="./family/data/pleiotropy/results_all_08212024.RData")


results_all = results_all[c(1,2,3,5,4,10,11,7),4:18]
results_all <- results_all %>% 
  mutate(EUR.pvalue = as.character(signif(EUR.pvalue,digits=3)),
         All.pvalue = as.character(signif(All.pvalue,digits=3))) 

results_all

results_all[,c(1,2,4:14)] = round(results_all[,c(1,2,4:14)],digits = 3)

write.csv(results_all,file="./family/data/pleiotropy/results_all_PGSmain_08212024.csv")



delta_all = delta_all[c(1,2,3,5,4,10,11,7),4:18]
delta_all <- delta_all %>% 
  mutate(EUR.pvalue = as.character(signif(EUR.pvalue,digits=3)),
         All.pvalue = as.character(signif(All.pvalue,digits=3))) 

delta_all

delta_all[,c(1,2,4:14)] = round(delta_all[,c(1,2,4:14)],digits = 3)

write.csv(delta_all,file="./family/data/pleiotropy/delta_all_nurture_08212024.csv")




#check the pearson correlations
load("./family/data/autism/autism_model_data_1000g_redo.RData")


trait=c("edu","schizophrenia","depression","bipolar","adhd","bmi","bmi_prive","bipolar1","bipolar2","neuroticism")
cor_pearson=data.frame(array(0,c(length(trait),4*4)))
colnames(cor_pearson)=c("cor","p","lwr","upr",
                        "cor_eur","p","lwr","upr","cor_afr","p","lwr","upr",
                        "cor_as","p","lwr","upr")
rownames(cor_pearson)=trait
h=0
for(i in trait){
  load(paste0("./family/data/pleiotropy/",i,"/",i,"_clean_trio_standardized.RData"))
  h=h+1
  data.model2=data.model1[-which(data.model1$race=="AMR"),]
  
  autism_prs=data.model$offspring_prs_normalized_1000g_combine[match(data.model2$offspring_id,data.model$offspring_id)]
  tmp=cor.test(autism_prs,data.model2$offspring_prs_standardized)
  cor_pearson[h,1:4]=c(tmp$estimate,tmp$p.value,tmp$conf.int[1],tmp$conf.int[2])
  
  
  data.model2=data.model1[which(data.model1$race=="EUR"),]
  autism_prs=data.model$offspring_prs_normalized_1000g_combine[match(data.model2$offspring_id,data.model$offspring_id)]
  tmp=cor.test(autism_prs,data.model2$offspring_prs_standardized)
  cor_pearson[h,5:8]=c(tmp$estimate,tmp$p.value,tmp$conf.int[1],tmp$conf.int[2])
  
  
  data.model2=data.model1[which(data.model1$race=="AFR"),]
  autism_prs=data.model$offspring_prs_normalized_1000g_combine[match(data.model2$offspring_id,data.model$offspring_id)]
  tmp=cor.test(autism_prs,data.model2$offspring_prs_standardized)
  cor_pearson[h,9:12]=c(tmp$estimate,tmp$p.value,tmp$conf.int[1],tmp$conf.int[2])
  
  data.model2=data.model1[which(data.model1$race=="SAS" | data.model1$race=="EAS"),]
  autism_prs=data.model$offspring_prs_normalized_1000g_combine[match(data.model2$offspring_id,data.model$offspring_id)]
  tmp=cor.test(autism_prs,data.model2$offspring_prs_standardized)
  cor_pearson[h,13:16]=c(tmp$estimate,tmp$p.value,tmp$conf.int[1],tmp$conf.int[2])

}
save(results_all,delta_all,cor_pearson,file="./family/data/pleiotropy/results_all.RData")


#Add AMR to the all population
#summary statistics can be found at PGS catalog
#or
#https://pgc.unc.edu/for-researchers/download-results/
setwd("./family/data/pleiotropy")
trait=c("edu","schizophrenia","depression","bipolar","adhd","bmi","bmi_prive","bipolar1","bipolar2","neuroticism")

results_all=array(0,c(length(trait),13))
results_all=data.frame(results_all)
rownames(results_all)=trait
colnames(results_all)=c("All.beta","All.se","All.pvalue",
                        "EUR.beta","EUR.se","EUR.pvalue",
                        "AFR.beta","AFR.se","AFR.pvalue",
                        "AS.beta","AS.se","AS.pvalue",
                        "#SNPs")
source("./family/R/PGScpt.R")
#N (All) = 1364; N (EUR) = 1250

delta_all=array(0,c(length(trait),12))
delta_all=data.frame(delta_all)
rownames(delta_all)=trait
#c("educational attainment","schizophrenia","depression","bipolar","ADHD","BMI","BMI_Prive")
colnames(delta_all)=c("All.beta","All.se","All.pvalue",
                      "EUR.beta","EUR.se","EUR.pvalue",
                      "AFR.beta","AFR.se","AFR.pvalue",
                      "AS.beta","AS.se","AS.pvalue")


h=0
for(i in trait){
  load(paste0("./family/data/pleiotropy/",i,"/",i,"_clean_trio_standardized.RData"))
  h=h+1
  data.model2=data.model1
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_nurture = T)
  results_all[h,1:3]=pgs_main$res_beta[c(1,2,4)]
  delta_all[h,1:3]=pgs_main$res_delta[c(1,2,4)]
  
  #EUR
  data.model2=data.model1[which(data.model1$race=="EUR"),]
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_nurture = T)
  results_all[h,4:6]=pgs_main$res_beta[c(1,2,4)]
  delta_all[h,4:6]=pgs_main$res_delta[c(1,2,4)]
  
  
  #AFR
  data.model2=data.model1[which(data.model1$race=="AFR"),]
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_nurture = T)
  results_all[h,7:9]=pgs_main$res_beta[c(1,2,4)]
  delta_all[h,7:9]=pgs_main$res_delta[c(1,2,4)]
  
  
  
  #AS
  data.model2=data.model1[which(data.model1$race=="SAS" | data.model1$race=="EAS"),]
  pgs_main=PGS.TRI(data.model2$offspring_prs_standardized,
                   data.model2$mother_prs_standardized,
                   data.model2$father_prs_standardized, 
                   GxE_int = FALSE,
                   parental_nurture = T)
  results_all[h,10:12]=pgs_main$res_beta[c(1,2,4)]
  delta_all[h,10:12]=pgs_main$res_delta[c(1,2,4)]
}

results_all$`#SNPs`=c(49265,585938,18746,741261,419357,1414,107018,  731828, 730189 ,53398)
results_all
save(results_all,delta_all,file="./family/data/pleiotropy/results_all_withAMR.RData")

