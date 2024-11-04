#Genotypic TDT test for rs2235370 and the other 8 SNPs in LD identified in GWAS Catalog of British and Chinese populations combined
#May 20, 2024
#----------------------------------------------------------------------------------

library(trio)
library(data.table)
setwd("./phewas/check_cis_snp/tdt/trio")



##############################################################################
##############################################################################
##############################################################################

#Prepare the data for genotypic TDT
target_snp = fread("OPGS007065.sscore.vars",header = F) #OMICSPRED is HG19 format
target_snp$location = gsub(".*:","",target_snp$V1)

#GWAS catalog is HG38 format, need to overlift to HG19
eur = fread("rs2235370_EUR03.tsv")
#liftover to hg19
eur_snp = eur$Location
eur_snp = paste0("chr",eur_snp)

eur_hg19 = fread("hglft_genome_EUR03.bed",header = F)
eur$Location = eur_hg19$V1
eur$loc = gsub(".*-","",eur$Location)
id = intersect(eur$loc,target_snp$location)
eur_intersect=eur[match(id,eur$loc),]


#Now check Asian
eas = fread("rs2235370_EAS.tsv")
#liftover to hg19
eas_snp = eas$Location
eas_snp = paste0("chr",eas_snp)

eas_hg19 = fread("hglft_genome_EAS.bed",header = F)
eas$Location = eas_hg19$V1
eas$loc = gsub(".*-","",eas$Location)
id = intersect(eas$loc,target_snp$location)
eas_intersect=eas[match(id,eas$loc),]


#all snps with r2 > 0.3
snp_rm = unique(c(eas_intersect$loc[which(eas_intersect$r2 > 0.3)],eur_intersect$loc)) #8 SNPs
snp_keep_name = target_snp$V1[match(snp_rm,target_snp$location)]


#save the 8 snps and rs2235370
#chr1:209946027-209946027 #hg19 position
snp_keep_name = c("1:209946027",snp_keep_name)
snp_keep_name = data.frame(snp_keep_name)
write.table(snp_keep_name,file="tdt/trio/snp_keep_name.txt",quote=F,row.names = F,col.names = F)



####################################################
#load raw file from PLINK
eas_trio = data.frame(fread("EAS_CLP.raw"))
eur_trio = data.frame(fread("EUR_CLP.raw"))
all_trio = data.frame(fread("ALL_CLP.raw"))


#load family information
#load("/dcs04/nilanjan/data/zwang/phewas/plink/oc/oc_trio.RData")
load("./phewas/check_cis_snp/tdt/trio/oc_trio.RData")
data.model1=data.model[complete.cases(data.model),]
table(data.model1$proband_race2)

#ASA  EUR 
#1126  778 

#First fit for all trios cross-population
data.model2=data.model[which(data.model$clp=="CL/P"),]
dim(data.model2)
#[1] 1466    8

colnames(all_trio)
res=list()
h=0
for(i in 7:15){
h=h+1
data_input_gtdt = data.frame(array(0,c(dim(data.model2)[1],3)))
data_input_gtdt[,1]=all_trio[match(data.model2$father_id,all_trio$IID),i]
data_input_gtdt[,2]=all_trio[match(data.model2$mother_id,all_trio$IID),i]
data_input_gtdt[,3]=all_trio[match(data.model2$offspring_id,all_trio$IID),i]
gtdt_input=as.vector(t(data_input_gtdt))
res[[h]] = tdt(gtdt_input,model="additive")
}
res_all = array(0,c(9,6))
for(i in 1:9){
  res_all[i,1] = res[[i]]$coef
  res_all[i,2] = res[[i]]$se
  res_all[i,3] = res[[i]]$RR
  res_all[i,4] = res[[i]]$lowerRR
  res_all[i,5] = res[[i]]$upperRR
  res_all[i,6] = res[[i]]$pval
}
colnames(res_all) = c("coef","se","RR","lowerRR","upperRR","P")
all_trio = fread("ALL_CLP.raw")
rownames(res_all) = colnames(all_trio)[7:15]

#EUR
data.model2=data.model1[which(data.model1$proband_race2=="EUR" & data.model1$clp=="CL/P"),]

res=list()
h=0
for(i in 7:15){
  h=h+1
  data_input_gtdt = data.frame(array(0,c(dim(data.model2)[1],3)))
  data_input_gtdt[,1]=eur_trio[match(data.model2$father_id,eur_trio$IID),i]
  data_input_gtdt[,2]=eur_trio[match(data.model2$mother_id,eur_trio$IID),i]
  data_input_gtdt[,3]=eur_trio[match(data.model2$offspring_id,eur_trio$IID),i]
  gtdt_input=as.vector(t(data_input_gtdt))
  res[[h]] = tdt(gtdt_input,model="additive")
}
res_eur = array(0,c(9,6))
for(i in 1:9){
  res_eur[i,1] = res[[i]]$coef
  res_eur[i,2] = res[[i]]$se
  res_eur[i,3] = res[[i]]$RR
  res_eur[i,4] = res[[i]]$lowerRR
  res_eur[i,5] = res[[i]]$upperRR
  res_eur[i,6] = res[[i]]$pval
}
colnames(res_eur) = c("coef","se","RR","lowerRR","upperRR","P")
rownames(res_eur) = rownames(res_all)



#Asian
data.model2=data.model1[which(data.model1$proband_race2=="ASA" & data.model1$clp=="CL/P"),]

res=list()
h=0
for(i in 7:15){
  h=h+1
  data_input_gtdt = data.frame(array(0,c(dim(data.model2)[1],3)))
  data_input_gtdt[,1]=eas_trio[match(data.model2$father_id,eas_trio$IID),i]
  data_input_gtdt[,2]=eas_trio[match(data.model2$mother_id,eas_trio$IID),i]
  data_input_gtdt[,3]=eas_trio[match(data.model2$offspring_id,eas_trio$IID),i]
  gtdt_input=as.vector(t(data_input_gtdt))
  res[[h]] = tdt(gtdt_input,model="additive")
}
res_eas = array(0,c(9,6))
for(i in 1:9){
  res_eas[i,1] = res[[i]]$coef
  res_eas[i,2] = res[[i]]$se
  res_eas[i,3] = res[[i]]$RR
  res_eas[i,4] = res[[i]]$lowerRR
  res_eas[i,5] = res[[i]]$upperRR
  res_eas[i,6] = res[[i]]$pval
}
colnames(res_eas) = c("coef","se","RR","lowerRR","upperRR","P")
rownames(res_eas) = rownames(res_all)

#save the results
save(res_eur,res_eas,res_all,file="gTDT_results.RData")




#To create results in Table 2
load("./phewas/check_cis_snp/tdt/trio/gTDT_results.RData")
order = c("1:209946027_A",
"1:210022901_T",
"1:209891906_A",
"1:209980757_A",
"1:210030755_G",
"1:209959614_T",
"1:209872606_T",
"1:209869121_C",
"1:209867116_A")

res_eur = data.frame(res_eur[match(order,rownames(res_eur)),c(3:6)])
res_eas = data.frame(res_eas[match(order,rownames(res_eas)),c(3:6)])
res_all = data.frame(res_all[match(order,rownames(res_all)),c(3:6)])

res_print = rbind(res_eur,res_eas,res_all)
res_print$population = c(rep("EUR",9),rep("Asian",9),rep("cross-pop",9))

write.csv(res_print,file="results_tableS16.csv")


maf_eas = fread("../tdt_EAS_CLP.FRQ")
maf_eur = fread("../tdt_EUR_CLP.FRQ")
maf_eas$name = paste0(maf_eas$SNP,"_",maf_eas$A1)
maf_eur$name = paste0(maf_eur$SNP,"_",maf_eur$A1)

maf_1 = maf_eur[match(rownames(res_eur),maf_eur$name),]
maf_2 = maf_eas[match(rownames(res_eas),maf_eas$name),]
maf_1$MAF_Asian = maf_2$MAF
maf_1 = maf_1[,-6]
colnames(maf_1)[5] = "MAF_EUR"
write.csv(maf_1,file="MAF_GENEVA.csv")


res_print[,1:3] = round(res_print[,1:3],digits = 2)
res_print$rr_ci = paste0("(",res_print$lowerRR,",",res_print$upperRR,")")
res_print[,4] = format(res_print[,4], digits=2)

res_print$population = c(rep("EUR",9),rep("Asian",9),rep("cross-pop",9))

write.csv(res_print,file="results_table2.csv")





################################################################################
#Re-generate the tables by standardizing the log RR using the genotype data SD
#To create results in Table 2
load("./phewas/check_cis_snp/tdt/trio/gTDT_results.RData")
order = c("1:209946027_A",
          "1:210022901_T",
          "1:209891906_A",
          "1:209980757_A",
          "1:210030755_G",
          "1:209959614_T",
          "1:209872606_T",
          "1:209869121_C",
          "1:209867116_A")

res_eur = data.frame(res_eur[match(order,rownames(res_eur)),])
res_eas = data.frame(res_eas[match(order,rownames(res_eas)),])
res_all = data.frame(res_all[match(order,rownames(res_all)),])


#standardize the RR and SD by genotype variance using MAF

maf_eas = fread("../tdt_EAS_CLP.FRQ")
maf_eur = fread("../tdt_EUR_CLP.FRQ")
maf_all = fread("../tdt_ALL_CLP.FRQ")
maf_eas$name = paste0(maf_eas$SNP,"_",maf_eas$A1)
maf_eur$name = paste0(maf_eur$SNP,"_",maf_eur$A1)
maf_all$name = paste0(maf_all$SNP,"_",maf_all$A1)



maf_1 = maf_eur[match(rownames(res_eur),maf_eur$name),]
maf_2 = maf_eas[match(rownames(res_eas),maf_eas$name),]
maf_3 = maf_all[match(rownames(res_all),maf_all$name),]


maf_1$MAF_Asian = maf_2$MAF
maf_1$MAF_all = maf_3$MAF
maf_1 = maf_1[,-6]
colnames(maf_1)[5] = "MAF_EUR"

res_eur$coef_standardized = res_eur$coef * sqrt(2*maf_1$MAF_EUR*(1-maf_1$MAF_EUR))
res_eur$se_standardized = res_eur$se * sqrt(2*maf_1$MAF_EUR*(1-maf_1$MAF_EUR))

res_eas$coef_standardized = res_eas$coef * sqrt(2*maf_1$MAF_Asian*(1-maf_1$MAF_Asian))
res_eas$se_standardized = res_eas$se * sqrt(2*maf_1$MAF_Asian*(1-maf_1$MAF_Asian))

res_all$coef_standardized = res_all$coef * sqrt(2*maf_1$MAF_all*(1-maf_1$MAF_all))
res_all$se_standardized = res_all$se * sqrt(2*maf_1$MAF_all*(1-maf_1$MAF_all))

res_eas$RR_standardized = exp(res_eas$coef_standardized)
res_eur$RR_standardized = exp(res_eur$coef_standardized)
res_all$RR_standardized = exp(res_all$coef_standardized)

res_eas$lowerRR_standardized = exp(res_eas$coef_standardized - qnorm(0.975)*res_eas$se_standardized)
res_eas$upperRR_standardized = exp(res_eas$coef_standardized + qnorm(0.975)*res_eas$se_standardized)

res_eur$lowerRR_standardized = exp(res_eur$coef_standardized - qnorm(0.975)*res_eur$se_standardized)
res_eur$upperRR_standardized = exp(res_eur$coef_standardized + qnorm(0.975)*res_eur$se_standardized)

res_all$lowerRR_standardized = exp(res_all$coef_standardized - qnorm(0.975)*res_all$se_standardized)
res_all$upperRR_standardized = exp(res_all$coef_standardized + qnorm(0.975)*res_all$se_standardized)



res_print = rbind(res_eur,res_eas,res_all)
res_print$population = c(rep("EUR",9),rep("Asian",9),rep("cross-pop",9))

res_print[,c(1:5,7:11)] = round(res_print[,c(1:5,7:11)],digits = 2)
res_print$rr_ci = paste0("(",res_print$lowerRR,",",res_print$upperRR,")")
res_print$rr_ci_standardized = paste0("(",res_print$lowerRR_standardized,",",res_print$upperRR_standardized,")")

res_print[,6] = format(res_print[,6], digits=2)

res_print$population = c(rep("EUR",9),rep("Asian",9),rep("cross-pop",9))

write.csv(res_print,file="results_table2_standardized.csv")
