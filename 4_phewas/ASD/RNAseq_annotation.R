#RNAseq results of the PGS main and indirect effect for ASD
#Manuscript results


library(biomaRt)
library(dplyr)
library(stringr)



#------------------------------------------------------------------
type="RNAseq"
load(paste0("./phewas/summary/redo/",type,"_results_smallsample.RData"))
res_pgs_main=data.frame(res_pgs_main)
res_indirect=data.frame(res_indirect)
info=fread(paste0("./phewas/data/","RNAseq_trait_validation_results_with_OMICSPRED_ID.csv"))
rm_name=unique(c(info$`OMICSPRED ID`[which(info$Internal_R2 < 0.1)],info$`OMICSPRED ID`[which(info$`#SNP` < 5)]))
#remove those with SNPs<5
res_pgs_main=res_pgs_main[-match(intersect(rm_name,rownames(res_pgs_main)),rownames(res_pgs_main)),]
res_indirect=res_indirect[-match(intersect(rm_name,rownames(res_indirect)),rownames(res_indirect)),]
res_pgs_main$ID=rownames(res_pgs_main)
res_indirect$ID=rownames(res_indirect)
res_info=merge(res_pgs_main,info,by.x = "ID",by.y="OMICSPRED ID")


## search by ensembl_gene_id

## first search by uniprot_gn_id in the ensembl database
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
searchDatasets(mart = ensembl, pattern = "hsapiens")
##                  dataset              description    version
## 80 hsapiens_gene_ensembl Human genes (GRCh38.p14) GRCh38.p14

allgene <- unique(info$`Ensembl ID`)
se.dat <- getBM(attributes = c("ensembl_gene_id","chromosome_name", "transcription_start_site"),
                filters = "ensembl_gene_id",
                values = allgene,
                mart = ensembl)
se.dat <- se.dat[se.dat$chromosome_name %in% c(as.character(1:22)), ]
se.dat <- se.dat %>% group_by(ensembl_gene_id) %>% summarise(chromosome_name=unique(chromosome_name), transcription_start_site=min(transcription_start_site))
#tmp <- allgene[!(allgene %in% se.dat$ensembl_gene_id)] # genes which have no results by ensembl_gene_id
#prot.anno <- inner_join(res_info[!(res_info$`Ensembl ID` %in% tmp),], se.dat, by=c("Ensembl ID"="ensembl_gene_id"))
prot.anno <- inner_join(res_info, se.dat, by=c("Ensembl ID"="ensembl_gene_id"))

#prot.anno.redo = merge(res_info,prot.anno[,c("ID","chromosome_name","transcription_start_site")],by="ID")
save(prot.anno,file="RNAseq_maineffect_anno_joint.RData")

indirect.anno=merge(res_indirect,prot.anno[,-c(2:19)],by="ID")
dim(indirect.anno)
#[1] 4966   29
save(indirect.anno,file="RNAseq_indirect_anno_joint.RData")

