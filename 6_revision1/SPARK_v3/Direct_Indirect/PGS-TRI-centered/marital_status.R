#extract UKB data for marital status

df.marital <- read.table(pipe("cut -f1,5848,5849,5850,5851,5852 /dcl01/arking/data/UK_Biobank/static/Phenotype/Downloads/tab_delim/ukb45830.tab"))
saveRDS(df.marital,file="/dcs04/nilanjan/data/zwang/family/pleiotropy/UKB/marital_ukbiobank.rds")


#test for male/female difference among married people 
marital <- readRDS("./family/data/pleiotropy/UKB/marital_ukbiobank.rds")
colnames(marital) = marital[1,]
marital = marital[-1,]
table(marital$f.6141.0.0)
table(marital$f.6141.0.1)
table(marital$f.6141.0.2)
table(marital$f.6141.0.3)
table(marital$f.6141.0.4)
#you can see that only f.6141.0.0 contains value 1 (living with Husband, wife or partner)
id_married = marital$f.eid[which(marital$f.6141.0.0 == "1")]


load("./family/data/pleiotropy/UKB/pgs_all_UKB_PCstandardized.RData") 
dat_UKB = dat_UKB[which(dat_UKB$IID %in% id_married),]
dim(dat_UKB)
trait=c("Educational Attainment","Schizophrenia","Depression","Bipolar Disorder","ADHD","BMI","Bipolar1","Bipolar2","Neuroticism","Insomnia","Chronotype","ASD")
table(dat_UKB$Sex)


res_ttest = list()
res_kstest = list()
for(i in 1:length(trait)){
  res_ttest[[i]] = t.test(dat_UKB[which(dat_UKB$Sex == "Female"),(2*i+27)],dat_UKB[which(dat_UKB$Sex == "Male"),(2*i+27)])
  
  res_kstest[[i]] = ks.test(dat_UKB[which(dat_UKB$Sex == "Female"),(2*i+27)],dat_UKB[which(dat_UKB$Sex == "Male"),(2*i+27)])$p
}

diff_var = function(x,y){
  diff = (mean(x)-mean(y))/sqrt((var(x)*(length(x)-1) + var(y)*(length(y)-1))/(length(x)+length(y)-2))
  sd = sqrt( (1/length(x) + 1/length(y))+ diff^2/(2*(length(x) + length(y))))
  res = list(diff = diff,sd =sd )
  return(res)
  
}

res_diff = list()
for(i in 1:length(trait)){
  res_diff[[i]] = diff_var(dat_UKB[which(dat_UKB$Sex == "Female"),(2*i+27)],dat_UKB[which(dat_UKB$Sex == "Male"),(2*i+27)])
}
res_table = do.call(rbind.data.frame, res_diff)
rownames(res_table) = trait
write.csv(res_table,file="./family/data/pleiotropy/UKB/diff_var_married.csv")
