#prepare the PLINK weight file to calculate gtex scores
#part of the code from https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R

# allele.qc = function(a1,a2,ref1,ref2) {
#   a1 = toupper(a1)
#   a2 = toupper(a2)
#   ref1 = toupper(ref1)
#   ref2 = toupper(ref2)
#   
#   ref = ref1
#   flip = ref
#   flip[ref == "A"] = "T"
#   flip[ref == "T"] = "A"
#   flip[ref == "G"] = "C"
#   flip[ref == "C"] = "G"
#   flip1 = flip
#   
#   ref = ref2
#   flip = ref
#   flip[ref == "A"] = "T"
#   flip[ref == "T"] = "A"
#   flip[ref == "G"] = "C"
#   flip[ref == "C"] = "G"
#   flip2 = flip;
#   
#   snp = list()
#   snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
#   snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
#   snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
#   snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
#   
#   return(snp)
# }

library(data.table)
#weight lists
names = fread("/dcs04/nilanjan/data/zwang/family/gtex/data/Tissue_list_GTex_V8_analyze.txt",header = F)
for(i in 1:dim(names)[1]){
wt.list = fread(paste0("/dcs04/nilanjan/data/zwang/family/gtex/data/weights/",names[i],".pos"))
#bim_SPARK = fread("/dcs05/ladd/NDEpi/data/Projects/InProgress/ziqiao/data/mergefinal.bim")
#load .RDat file
for(j in 1:dim(wt.list)[1]){
load(paste0("/dcs04/nilanjan/data/zwang/family/gtex/data/weights/",wt.list$WGT[j]))
#find the best set of weights
best.mod = which.max(cv.performance["rsq",])
info_score = cbind(snps[,c("V2","V5")],wgt.matrix[match(snps$V2,rownames(wgt.matrix)),best.mod])
#save this file for calculating scores in plink
write.table(info_score,file=paste0("/dcs04/nilanjan/data/zwang/family/gtex/data/weights_clean/",names[i],"/",wt.list$ID[j],".txt"),row.names = F,col.names = F,quote = F)
}
}