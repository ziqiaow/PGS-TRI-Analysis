#Summarize the results for GxE analysis
#supplementary table

setwd("~/work/family/data/autism/v3/GxE/homo")
load("interaction_results_clean.RData")
name = rownames(tmp_interaction_all)
order = c(1,6:11,2,5,12,13,14,15,18:22,3,23,24,25,16,17)

tmp_interaction_all = tmp_interaction_all[name[order],]
tmp_interaction_EUR = tmp_interaction_EUR[name[order],]
tmp_interaction_AMR = tmp_interaction_AMR[name[order],]

tmp_interaction_AFR_match = tmp_interaction_AFR[match(name[order],rownames(tmp_interaction_AFR)),]
rownames(tmp_interaction_SAS)[14] = "PGS x factor(smoke_tob_ever_preg)1"
tmp_interaction_SAS_match = tmp_interaction_SAS[match(name[order],rownames(tmp_interaction_SAS)),]

tmp_int = cbind(tmp_interaction_all,tmp_interaction_EUR,tmp_interaction_AFR_match,tmp_interaction_AMR,tmp_interaction_SAS_match)
colnames(tmp_int) = paste0(colnames(tmp_interaction_all),rep(c("All","EUR","AFR","AMR","EAS/SAS"),each = 3))
write.csv(tmp_int,file="gxe_homo_all_paper.csv")
