
#m.out.mahal.exact and m.out.mahal performed the best
load("./family/simulation/UKB/match_ps_mahal_exact_geoonly.RData")
load("./family/simulation/UKB/selected_final.RData")
#let's check the pairs
id_match = m.out.mahal.exact$match.matrix
id_parental = data.frame(array(0,c(dim(id_match)[1],2)))
id_parental[,1] = rownames(id_match)
id_parental[,2] = id_match[,1]
id_parental = id_parental[complete.cases(id_parental),]
dim(id_parental)
#[1] 150253      2
selected_match = data.frame(selected_final[,c("ID","birth.x","birth.y","assessment_center","education_match","sex","education","townsend_index","BMI","age","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")])
selected_match$sex = as.numeric(selected_match$sex)
rownames(selected_match) = selected_match$ID
table(selected_match$education)
class(selected_match$education)
selected_match$education = as.numeric(selected_match$education)
selected_match$education[selected_match$education<0] = 0
pheno_father = selected_match[match(id_parental[,1],selected_match$ID),]
pheno_mother = selected_match[match(id_parental[,2],selected_match$ID),]
save(id_parental,pheno_father,pheno_mother,selected_match,file="match_geo/matched_parents_pheno.RData")
id_total = c(id_parental[,1],id_parental[,2])
id_total = data.frame(cbind(id_total,id_total))
colnames(id_total) = c("#FID","IID")
write.table(id_total,file="match_geo/plink/keep_parents.txt",row.names = F,col.names = T,quote = F)

