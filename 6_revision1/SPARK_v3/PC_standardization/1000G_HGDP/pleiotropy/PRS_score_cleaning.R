#----------------------------------
#genotype data cleaning
#GEARs
#March 29, 2023
#----------------------------------
#https://choishingwan.github.io/PRS-Tutorial/base/
#Complete steps:
#1. Download PGS file from PGS catalog
#2. Remove ambiguous snps and missing harmonized position for HG38 (note, use hm_pos and hm_chr in the file for harmonized file)
#3. Save this new snp file and extract them from the complete data using plink, the command is in C:/Users/ziqia/Desktop/work/family/data/autism/plink/extract_snp.sh
#4. QC the extracted files, the plink command is C:/Users/ziqia/Desktop/work/family/data/autism/plink/qc_snp.sh
#5. Merge the extracted bim files, the plink command is C:/Users/ziqia/Desktop/work/family/data/autism/plink/combine_bim2.sh
#6. Read in the combined bim file and identify SNPs that require strand flipping, save the score file of effect allele, rsid, and weight
#7. Use plink to calculate the PRS score, the plink command is C:/Users/ziqia/Desktop/work/family/data/autism/plink/PRS_plink_merged.sh
#Note that the plink actually identify the number of effect alleles in each individual, therefore, no need to change the weight or major/minor allele information
#however, if you insist on doing so, you can find these steps at the end of this R file. I've calculated using plink2 and plink1.9(better documentation of how they calculate PRS), using both matched and unmatched bim/bed/fam files, results are perfectly the same (cor=1).

setwd("C:/Users/ziqia/Desktop/work/family/data/pleiotropy")
#load bim file for SPARK and 1KG (the reason using overlapping SNPs for PC-standardization and 1KG mean variance standardization, so PRS are based on the same SNPs)
load("C:/Users/ziqia/Desktop/work/family/data/autism/v3/plink/bim_all_match1KG.RData")
bimfile_spark = data.frame(bimfile_spark)
#first get the SNPs from PGS catalog
library(data.table)

trait=c("edu","schizophrenia","depression","bipolar","bmi_prive","bipolar1","bipolar2","neuroticism","insomnia","chronotype")
path0 <- "C:/Users/ziqia/Desktop/work/family/data/pleiotropy/"
path1 <- "C:/Users/ziqia/Desktop/work/family/data/pleiotropy/redo/" #for ADHD

snp=list()
h=0
files = list()
for(i in trait){
  h=h+1
  
  files[[h]] <- list.files(path = paste0(path0,i), pattern = "*hmPOS_GRCh38.txt$")
  snp[[h]] = fread(paste0(path0,i,"/",files[[h]]),header = T)
  
}
for(i in 1:length(trait)){
  print(i)
  print(colnames(snp[[i]]))
  if( !any((colnames(snp[[i]]) %in% "other_allele") == "TRUE")){
    snp[[i]]$other_allele = snp[[i]]$hm_inferOtherAllele
  }
  print(i)
  print(colnames(snp[[i]]))
}
for(i in 1:length(trait)){
  tmp = snp[[i]]
  print(i)
  print(dim(tmp))
  rm_id=which( (tmp$effect_allele=="A" & tmp$other_allele=="T") |
                 (tmp$effect_allele=="T" & tmp$other_allele=="A") |
                 (tmp$effect_allele=="C" & tmp$other_allele=="G") |
                 (tmp$effect_allele=="G" & tmp$other_allele=="C") )
  if(length(rm_id) > 1){
    snp[[i]]=tmp[-rm_id,]
    
  }
  print(dim(snp[[i]]))
 
  
  
  
}


# for(i in trait){
#   h=h+1
#   if(i != "ADHD"){
#     files <- list.files(path = paste0(path0,i), pattern = "*hmPOS_GRCh38.txt$")
#     snp[[h]] = fread(paste0(path0,i,"/",files),header = T)
#   } else {
#     files <- list.files(path = paste0(path1,i), pattern = "*hmPOS_GRCh38.txt$")
#     snp[[h]] = fread(paste0(path1,i,"/",files),header = T)}
#   
# }
#   


#then we combine the BIM files with score file

#match with our PRS weight file
for(i in 1:length(trait)){
tmp = data.frame(snp[[i]])
score = merge(bimfile_spark,tmp,by.x = c("chr","position"),by.y = c("hm_chr","hm_pos"))
score_output = score[,c("rsid","effect_allele","effect_weight")]
score_output = unique(score_output)
print(dim(score_output))
#name = gsub("_.*$", "", files[[i]] )
write.table(score_output,file=paste0("C:/Users/ziqia/Desktop/work/family/data/pleiotropy/redo/plink/",trait[i],"_score.txt"),row.names = F,col.names = F,quote = F)



}

# [1] 48546     3
# [1] 574170      3
# [1] 18253     3
# [1] 914179      3
# [1] 105350      3
# [1] 903999      3
# [1] 902155      3
# [1] 52687     3
# [1] 886084      3
# [1] 913313      3

#type this in linux to create pgs_trait.txt file
#ls > /dcs04/nilanjan/data/zwang/family/pleiotropy/redo/pgs_trait.txt