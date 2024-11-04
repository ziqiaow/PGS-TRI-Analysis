#Simulations for family-based study
#Use propensity score matching to find between male and female partners in the UK Biobank
#04/21/2024



#-----------------------------------------------------------
setwd("./family/simulation/UKB")
library(MatchIt)
library(data.table)
#we first load the UKB data
selected_final = readRDS('./CCA/data/t2d/UKB/pheno_unrelated_clean_t2d_03032024.rds')
reclean = readRDS("./CCA/data/UKB/rotation/T2D_CAD_clean_04162024.rds")
good_ids <- fread("./CCA/data/UKB/british_unrelated_good_qc_ids")
setnames(good_ids,"ID")
covariates <- fread("./CCA/data/t2d/UKB/covariates", na.strings = c("",".","NA", "Prefer not to answer","Do not know"))
covariates=covariates[,-c(20,21)]
merged <- Reduce(function(x,y) merge(x,y, by = "ID"), list(good_ids,reclean,covariates))

merged[,assessment_center:=factor(assessment_center)]

selected_final = selected_final[match(merged$ID,selected_final$ID),]
merged$birth.x = selected_final$birth.x
merged$birth.y = selected_final$birth.y

selected_final=merged


#use sex as the treatment group for matching, covariates including birth.x, birth.y, assessment centers (exact match), educational attainment for matching
selected_final$sex[selected_final$sex == "Male"] = 0
selected_final$sex[selected_final$sex == "Female"] = 1
table(selected_final$sex)

#0      1 
#156203 181183 

#Educational attainment: in UKB, field 6138 Qualifications
#https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=6138
#check when I did the GxE paper for breast cancer social economic phenotypes

#6138 Qualifications (education level)
#1	College or University degree
# 2	A levels/AS levels or equivalent
# 3	O levels/GCSEs or equivalent
# 4	CSEs or equivalent
# 5	NVQ or HND or HNC or equivalent
# 6	Other professional qualifications eg: nursing, teaching
# -7	None of the above
# -3	Prefer not to answer

social_ukbiobank=readRDS("./GXE/ukbiobank/data/social_ukbiobank.rds")
colnames(social_ukbiobank)=social_ukbiobank[1,]
social_ukbiobank=social_ukbiobank[-1,]
social_ukbiobank=social_ukbiobank[,-c(12:26)]
colnames(social_ukbiobank)=c("IID","household_income","work_time","age.education","moodswing","education","employment","employment_correct","depression","neuroticism","mental_health")

selected_final$education = social_ukbiobank$education[match(as.character(selected_final$ID),as.character(social_ukbiobank$IID))]
selected_final$education_match = NA
selected_final$education_match[which(selected_final$education==1)]="high"
selected_final$education_match[which(selected_final$education== -7 | selected_final$education== 4 |selected_final$education== 5 |selected_final$education== 6)]="low"
selected_final$education_match[which(selected_final$education== 2 | selected_final$education== 3  )]="intermediate"
selected_final$education_match=factor(selected_final$education_match,levels=c("low","intermediate","high"))

##22189 Townsend deprivation index at recruitment
#check this html file for Townsend deprivation index code in the data
#/dcs05/legacy-dcl01-arking/data/UK_Biobank/static/Phenotype/Documentation/ukb45830.html
#column 457, code 189-0.0 (different from UKB website)
# header.ph=readRDS("./GXE/ukbiobank/data/pheno_header.rds")
# col=c(grep("f.eid",header.ph),grep("f.189.0.0",header.ph))
# paste0(col,collapse=",")
# #1,457
# #3. read from server with a small memory size is fine
# df.social <- read.table(pipe("cut -f1,457 /dcl01/arking/data/UK_Biobank/static/Phenotype/Downloads/tab_delim/ukb45830.tab"))
# saveRDS(df.social,file="/users/zwang4/ukbiobank/townsend_index_ukb.rds")

townsend <- readRDS("plink/townsend_index_ukb.rds")
colnames(townsend)=townsend[1,]
townsend=townsend[-1,]
colnames(townsend)=c("IID","townsend_index")
townsend$townsend_index = as.numeric(townsend$townsend_index)
selected_final$townsend_index = townsend$townsend_index[match(as.character(selected_final$ID),as.character(townsend$IID))]
save(selected_final,file="selected_final.RData")


library(ggfortify)
selected_final_plot = selected_final[,c("birth.x","birth.y","education","education_match","townsend_index")]
selected_final_plot[selected_final_plot$education == -3] = NA
selected_final_plot = selected_final_plot[complete.cases(selected_final_plot),]

autoplot(kmeans(selected_final_plot[,c("birth.x","birth.y")], 100), data = selected_final_plot[,c("birth.x","birth.y")], frame = TRUE)
pc_ea = prcomp(selected_final_plot[,c("birth.x","birth.y")],scale. = T)

autoplot(pc_ea, data = selected_final_plot, colour = 'education_match')
autoplot(pc_ea, data = selected_final_plot, colour = 'education')
autoplot(pc_ea, data = selected_final_plot, colour = 'townsend_index')

require(ggplot2)
qplot(birth.x, birth.y, data = selected_final_plot, colour = factor(education))
qplot(birth.x, birth.y, data = selected_final_plot, colour = townsend_index)
qplot(birth.x, birth.y, data = selected_final_plot, colour = education_match)



selected_match = data.frame(selected_final[,c("ID","birth.x","birth.y","assessment_center","education_match","sex","education")])
selected_match = selected_match[complete.cases(selected_match),]
class(selected_match$sex)
selected_match$sex = as.numeric(selected_match$sex)
rownames(selected_match) = selected_match$ID
table(selected_match$education)
class(selected_match$education)
selected_match$education = as.numeric(selected_match$education)
selected_match$education[selected_match$education<0] = 0

match.out <- matchit(sex ~ birth.x + birth.y + education_match, 
                     data = selected_match, method = "nearest", exact = ~assessment_center) 
summary(match.out)
# Sample Sizes:
#   Control Treated
# All        151321  175206
# Matched    150946  150946
# Unmatched     375   24260
# Discarded       0       0

match.out2 <- matchit(sex ~ birth.x + birth.y + education, 
                     data = selected_match, method = "nearest", distance = "glm", exact = ~assessment_center) 
summary(match.out2)
save(match.out2,file="match_ps.RData")
plot(match.out2, type = "jitter", interactive = FALSE)  

# Mahalanobis matching within propensity score calipers
m.out.mahal <- matchit(sex ~ birth.x + birth.y + education + assessment_center, 
                       data = selected_match, method = "nearest", 
                       caliper = .1, mahvars = ~ birth.x + birth.y + education) 
summary(m.out.mahal)
save(m.out.mahal,file="match_ps_mahal.RData")



m.out.mahal.noas <- matchit(sex ~ birth.x + birth.y + education, 
                       data = selected_match, method = "nearest", 
                       caliper = .1, mahvars = ~ birth.x + birth.y + education) 
summary(m.out.mahal.noas)
save(m.out.mahal.noas,file="match_ps_mahal_noas.RData")


# Mahalanobis matching within propensity score calipers
m.out.mahal.exact <- matchit(sex ~ birth.x + birth.y + education, 
                       data = selected_match, method = "nearest", 
                       caliper = .1, mahvars = ~ birth.x + birth.y + education, exact = ~assessment_center) 
summary(m.out.mahal.exact)
save(m.out.mahal.exact,file="match_ps_mahal_exact.RData")



# Mahalanobis matching within propensity score calipers, geographical regions only
m.out.mahal.exact <- matchit(sex ~ birth.x + birth.y , 
                             data = selected_match, method = "nearest", 
                             caliper = .1, mahvars = ~ birth.x + birth.y, exact = ~assessment_center) 
summary(m.out.mahal.exact)
save(m.out.mahal.exact,file="match_ps_mahal_exact_geoonly.RData")





#Summarize balance statistics in a plot
plot(summary(m.out.mahal.noas),xlim=c(0,0.4))

#Summarize balance statistics in a plot
plot(summary(m.out.mahal),xlim=c(0,0.4))

#Summarize balance statistics in a plot
plot(summary(match.out2),xlim=c(0,0.4))

#Summarize balance statistics in a plot
plot(summary(m.out.mahal.exact),xlim=c(0,0.3))


#m.out.mahal.exact and m.out.mahal performed the best
load("./family/simulation/UKB/match_ps_mahal_exact.RData")
#let's check the pairs
id_match = m.out.mahal.exact$match.matrix
id_parental = data.frame(array(0,c(dim(id_match)[1],2)))
id_parental[,1] = rownames(id_match)
id_parental[,2] = id_match[,1]
id_parental = id_parental[complete.cases(id_parental),]
dim(id_parental)
#[1] 145956      2
selected_match = data.frame(selected_final[,c("ID","birth.x","birth.y","assessment_center","education_match","sex","education","townsend_index","BMI","age","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")])
selected_match$sex = as.numeric(selected_match$sex)
rownames(selected_match) = selected_match$ID
table(selected_match$education)
class(selected_match$education)
selected_match$education = as.numeric(selected_match$education)
selected_match$education[selected_match$education<0] = 0
pheno_father = selected_match[match(id_parental[,1],selected_match$ID),]
pheno_mother = selected_match[match(id_parental[,2],selected_match$ID),]
save(id_parental,pheno_father,pheno_mother,selected_match,file="matched_parents_pheno.RData")
id_total = c(id_parental[,1],id_parental[,2])
id_total = data.frame(cbind(id_total,id_total))
colnames(id_total) = c("#FID","IID")
write.table(id_total,file="plink/keep_parents.txt",row.names = F,col.names = T,quote = F)

