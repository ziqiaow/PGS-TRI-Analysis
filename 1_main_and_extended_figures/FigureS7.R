#----------------------------------------------------------
#Extended Figure 7
#----------------------------------------------------------
library(patchwork)
library(dplyr, warn = FALSE)
library(ggplot2)


library(data.table)
setwd("./phewas/plink/oc/summary/redo_06122024")

#qqplot
qqplot0 <- function(pval, title, col_point="blue", cex_point = 7/12) {
  match.loose <- function(a, b) sapply(a, function(x) which.min(abs(x-b)))
  pval <- sort(pval[!is.na(pval)])
  n0 <- length(pval)
  funnelx <- seq(log10(1+n0), 0, by=-0.01)
  Select0 <- unique(match.loose(10^(-funnelx)*(1+n0), 1:n0))
  lobs <- -log10(pval[Select0])
  lexp <- -log10(Select0/(1+n0))
  funnell<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.975,x,n0+1-x)))
  funnelu<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.025,x,n0+1-x)))
  # Set minimal margins (bottom, left, top, right)
  par(mar=c(4, 3, 2, 0.5),lwd=0.1, tcl=-0.25, mgp=c(1.5, 0.5, 0))
  plot(lexp, lobs, xlab=expression(Expected ~ ~-log[10](italic(p))), ylab=expression(Observed ~ ~-log[10](italic(p))), main=title, pch=16, cex=7/12, col="red", xlim=c(0, max(lexp)), ylim=c(0,max(funnelu)), type="n", yaxt="n",xaxt="n", cex.main=7/12, cex.lab=7/12, cex.axis=7/12)
  axis(1, at = seq(0, 3,1), labels = seq(0, 3, 1), lwd = 0, lwd.ticks=0.1,cex.axis=7/12)
  axis(2, at = seq(0, 5,1), labels = seq(0, 5, 1), lwd = 0, lwd.ticks=0.1,cex.axis=7/12)
   polygon(c(funnelx,rev(funnelx)),c(funnell,rev(funnelu)),col="grey",border="grey")
  lines(c(0,max(lexp)),c(0,max(lexp)))
  points(lexp, lobs, pch=16, cex=cex_point, col=col_point)
  text(max(lexp)*0.7, 0, as.expression(substitute(lambda[GC]==ll, list(ll=round(qchisq(median(pval), 1, lower=F)/qchisq(0.5, 1), 3)))), pos=3, cex=7/12)

}



qqplot <- function(pval, title, col_point="blue", cex_point = 7/12) {
  match.loose <- function(a, b) sapply(a, function(x) which.min(abs(x-b)))
  pval <- sort(pval[!is.na(pval)])
  n0 <- length(pval)
  funnelx <- seq(log10(1+n0), 0, by=-0.01)
  Select0 <- unique(match.loose(10^(-funnelx)*(1+n0), 1:n0))
  lobs <- -log10(pval[Select0])
  lexp <- -log10(Select0/(1+n0))
  funnell<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.975,x,n0+1-x)))
  funnelu<-sapply(10^(-funnelx)*(1+n0), function(x) -log10(qbeta(0.025,x,n0+1-x)))
  
  # Set minimal margins (bottom, left, top, right)
  par(mar=c(4, 0.5, 2, 0.5),lwd=0.1, tcl=-0.25, mgp=c(1.5, 0.5, 0))
  
  plot(lexp, lobs, xlab=expression(Expected ~ ~-log[10](italic(p))), ylab="", main=title, 
       pch=16, cex=7/12, col="red", xlim=c(0, max(lexp)), ylim=c(0,max(funnelu)), 
       type="n", cex.main=7/12, cex.lab=7/12, cex.axis=7/12, yaxt="n",xaxt="n")
  axis(1, at = seq(0, 3,1), labels = seq(0, 3, 1), lwd = 0, lwd.ticks=0.1,cex.axis=7/12) #https://www.r-bloggers.com/2020/08/custom-tick-marks-with-rs-base-graphics-system/
  polygon(c(funnelx,rev(funnelx)),c(funnell,rev(funnelu)),col="grey",border="grey")
  lines(c(0,max(lexp)),c(0,max(lexp)))
  points(lexp, lobs, pch=16, cex=cex_point, col=col_point)
  text(max(lexp)*0.7, 0, as.expression(substitute(lambda[GC]==ll, list(ll=round(qchisq(median(pval), 1, lower=F)/qchisq(0.5, 1), 3)))), pos=3, cex=7/12)
  
}

# qqplot for OC

#------------------------------------------------------------------
type="RNAseq"
load("./phewas/plink/oc/summary/redo_06122024/RNAseq_results.RData")
res_nurture=data.frame(res_nurture)


load("./phewas/plink/oc/summary/pgsmainonly/RNAseq_results.RData")
res_pgs_main=data.frame(res_pgs_main)

info=fread(paste0("./phewas/data/","RNAseq_trait_validation_results_with_OMICSPRED_ID.csv"))
rm_name=unique(info$`OMICSPRED ID`[which(info$Internal_R2<0.1 | info$`#SNP` < 5)])
#remove those with SNPs<5
res_pgs_main=res_pgs_main[-match(intersect(rm_name,rownames(res_pgs_main)),rownames(res_pgs_main)),]
res_nurture=res_nurture[-match(intersect(rm_name,rownames(res_nurture)),rownames(res_nurture)),]


cairo_pdf("./family/manuscript/draft4/plots/extended/Extended_Figure7.pdf",width=180/25.4, height=120/25.4)
layout(matrix(1:6, 2, 3, byrow=TRUE), widths=c(1.2, 1, 1))
qqplot0(res_pgs_main$p_clp_all,title="All Population (#Trios = 1466); DE",col_point = "royal blue")
qqplot(res_pgs_main$p_clp_EUR,title="EUR (#Trios = 575); DE",col_point = "royal blue")
qqplot(res_pgs_main$p_clp_AS,title="Asian (#Trios = 891); DE",col_point = "royal blue")
qqplot0(res_nurture$p_clp_all,title="All Population (#Trios = 1466); \u03B4-IDE",col_point = "dark green")
qqplot(res_nurture$p_clp_EUR,title="EUR (#Trios = 575); \u03B4-IDE",col_point = "dark green")
qqplot(res_nurture$p_clp_AS,title="Asian (#Trios = 891); \u03B4-IDE",col_point = "dark green")
dev.off()



#test for heterogeneity for TRAF3IP3 between EUR and Asian
#paried comparisons
heterogeneity.pair=function(beta1,beta2,se1,se2){
  test.stat=(beta1-beta2)/sqrt(se1^2+se2^2)
  p=2 * pnorm(abs(test.stat), lower.tail = FALSE)
  print(paste0("heterogeneity test p-value is ",p))
}

heterogeneity.pair(res_pgs_main[3503,]$beta_clp_EUR,res_pgs_main[3503,]$beta_clp_AS,res_pgs_main[3503,]$se_clp_EUR,res_pgs_main[3503,]$se_clp_AS)
#[1] "heterogeneity test p-value is 0.0134867598819943"