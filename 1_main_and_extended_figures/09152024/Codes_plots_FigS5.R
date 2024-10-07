## Extended Fig 5

library(patchwork)
library(dplyr, warn = FALSE)
library(ggplot2)
library(data.table)

#qqplot
qqplot <- function(pval, title, col_point="blue", cex_point = 0.6) {
  match.loose <- function(a, b) sapply(a, function(x) which.min(abs(x-b)))
  pval <- sort(pval[!is.na(pval)])
  n0 <- length(pval)
  funnelx <- seq(log10(1+n0), 0, by=-0.01)
  Select0 <- unique(match.loose(10^(-funnelx)*(1+n0), 1:n0))
  lobs <- -log10(pval[Select0])
  lexp <- -log10(Select0/(1+n0))
  funnell<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.975,x,n0+1-x)))
  funnelu<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.025,x,n0+1-x)))
  plot(lexp, lobs, xlab=expression(Expected ~ ~-log[10](italic(p))), ylab=expression(Observed ~ ~-log[10](italic(p))), main=title, pch=16, cex=0.5, col="red", xlim=c(0, max(lexp)), ylim=c(0,max(funnelu)), type="n", cex.main=1, cex.lab=1, cex.axis=1.5)
  polygon(c(funnelx,rev(funnelx)),c(funnell,rev(funnelu)),col="grey",border="grey")
  lines(c(0,max(lexp)),c(0,max(lexp)))
  points(lexp, lobs, pch=16, cex=cex_point, col=col_point)
  text(max(lexp)*0.7, 0, as.expression(substitute(lambda[GC]==ll, list(ll=round(qchisq(median(pval), 1, lower=F)/qchisq(0.5, 1), 3)))), pos=3)
  
}



# qqplot for ASD


#------------------------------------------------------------------
type="RNAseq"
load(paste0("./phewas/summary/redo_06122024/",type,"_results.RData"))
res_nurture=data.frame(res_nurture)

load(paste0("./phewas/summary/pgsmainonly/redo_06122024/",type,"_results.RData"))
res_pgs_main=data.frame(res_pgs_main)

info=fread(paste0("./phewas/data/","RNAseq_trait_validation_results_with_OMICSPRED_ID.csv"))
rm_name=unique(info$`OMICSPRED ID`[which(info$Internal_R2 < 0.1 | info$`#SNP` < 5)])
#remove those with SNPs<5
res_pgs_main=res_pgs_main[-match(intersect(rm_name,rownames(res_pgs_main)),rownames(res_pgs_main)),]
res_nurture=res_nurture[-match(intersect(rm_name,rownames(res_nurture)),rownames(res_nurture)),]


pdf("./family/manuscript/plots/final/Extended_Figure5_09152024.pdf",width = 12)
par(mfrow=c(2,2))
qqplot(res_pgs_main$p_all,title="All (#Trios = 1517); DE",col_point = "royal blue")
qqplot(res_pgs_main$p_EUR,title="EUR (#Trios = 1250); DE",col_point = "royal blue")
qqplot(res_nurture$p_all,title="All (#Trios = 1517); \u03B4-IDE",col_point = "dark green")
qqplot(res_nurture$p_EUR,title="EUR (#Trios = 1250); \u03B4-IDE",col_point = "dark green")
dev.off()



png("./family/manuscript/plots/final/Extended_Figure5_09152024.png",width=8, height=6, units="in", res=320)
par(mfrow=c(2,2))
qqplot(res_pgs_main$p_all,title="All (#Trios = 1517); DE",col_point = "royal blue")
qqplot(res_pgs_main$p_EUR,title="EUR (#Trios = 1250); DE",col_point = "royal blue")
qqplot(res_nurture$p_all,title="All (#Trios = 1517); \u03B4-IDE",col_point = "dark green")
qqplot(res_nurture$p_EUR,title="EUR (#Trios = 1250); \u03B4-IDE",col_point = "dark green")
dev.off()

