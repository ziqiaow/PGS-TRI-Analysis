#03/01/2025
#Original data cleaning in C:\Users\ziqia\Desktop\work\family\data\pleiotropy\adhd_redo

#Prepare data for input in PLINK
#Combine the scoring files from 22 chromosomes to 1 file, save to snpid, effect allele, effect size
library(data.table)
score_file = fread("/dcs04/nilanjan/data/zwang/family/pleiotropy/redo/adhd/prscs_results/eur_pst_eff_a1_b0.5_phiauto_chr1.txt",header = F)
score_file = data.frame(score_file)
score_file = score_file[,c(2,4,6)]
for(i in 2:22){
  tmp = fread(paste0("/dcs04/nilanjan/data/zwang/family/pleiotropy/redo/adhd/prscs_results/eur_pst_eff_a1_b0.5_phiauto_chr",i,".txt"),header = F)
  tmp = data.frame(tmp)
  tmp = tmp[,c(2,4,6)]
  score_file = rbind(score_file,tmp)
}
write.table(score_file,file="/dcs04/nilanjan/data/zwang/family/pleiotropy/redo/adhd_score.txt",col.names = F,row.names = F,quote = F)

dim(score_file)
#[1] 1025311       3
