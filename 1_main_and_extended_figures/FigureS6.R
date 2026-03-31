#----------------------------------------------------------
#Extended Figure 6
#----------------------------------------------------------
library(patchwork)
library(dplyr, warn = FALSE)
library(ggplot2)
library(data.table)

#qqplot

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

#sample-size adjusted genomic inflation factor
qqplot_adjusted_lambda0 <- function(pval, title, N_case,N_control,col_point="blue", cex_point = 7/12) {
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
  
  plot(lexp, lobs, xlab=expression(Expected ~ ~-log[10](italic(p))), ylab=expression(Observed ~ ~-log[10](italic(p))), main=title, pch=16, cex=7/12, col="red", xlim=c(0, max(lexp)), ylim=c(0,max(funnelu)), type="n", cex.main=7/12, cex.lab=7/12, cex.axis=7/12,yaxt="n",xaxt="n")
  axis(1, at = seq(0, 3,1), labels = seq(0, 3, 1), lwd = 0, lwd.ticks=0.1,cex.axis=7/12)
  axis(2, at = seq(0, 5,1), labels = seq(0, 5, 1), lwd = 0, lwd.ticks=0.1,cex.axis=7/12)
  polygon(c(funnelx,rev(funnelx)),c(funnell,rev(funnelu)),col="grey",border="grey")
  lines(c(0,max(lexp)),c(0,max(lexp)))
  points(lexp, lobs, pch=16, cex=cex_point, col=col_point)
  text(max(lexp)*0.7, 0, as.expression(substitute(lambda[GC_adjusted]==ll, list(ll= round(1+500 * (qchisq(median(pval), 1, lower=F)/qchisq(0.5, 1) -1)*(1/N_case + 1/N_control) , 3)))), pos=3, cex=7/12)
  
}

qqplot_adjusted_lambda <- function(pval, title, N_case, N_control, col_point="blue", cex_point = 7/12) {
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
  axis(1, at = seq(0, 3,1), labels = seq(0, 3, 1), lwd = 0, lwd.ticks=0.1,cex.axis=7/12)
  polygon(c(funnelx,rev(funnelx)),c(funnell,rev(funnelu)),col="grey",border="grey")
  lines(c(0,max(lexp)),c(0,max(lexp)))
  points(lexp, lobs, pch=16, cex=cex_point, col=col_point)
  text(max(lexp)*0.7, 0, as.expression(substitute(lambda[GC_adjusted]==ll, list(ll= round(1+500 * (qchisq(median(pval), 1, lower=F)/qchisq(0.5, 1) -1)*(1/N_case + 1/N_control) , 3)))), pos=3, cex=7/12)
}

# qqplot for ASD


#------------------------------------------------------------------
type="RNAseq"
load(paste0("./family/data/phewas/results/",type,"_final_results_homo.RData"))
direct_effect=data.frame(direct_effect)
indirect_effect=data.frame(indirect_effect)


#load true SNP counts
load(paste0("./family/data/phewas/results/SNP_count/",type,"_final_SNPcount.RData"))
direct_effect$'#SNP_final' = count[match(rownames(direct_effect),names(count))]
indirect_effect$'#SNP_final' = count[match(rownames(indirect_effect),names(count))]


#remove those with #SNPs<5
rm_name1 = rownames(direct_effect)[which(direct_effect$`#SNP_final` < 5)]
info=fread(paste0("./phewas/data/",type,"_trait_validation_results_with_OMICSPRED_ID.csv"))
rm_name=info$`OMICSPRED ID`[which(info$Internal_R2<0.1 | info$`#SNP` < 5 )]
rm_name = unique(c(rm_name,rm_name1))

direct_effect=direct_effect[-match(intersect(rm_name,rownames(direct_effect)),rownames(direct_effect)),]
indirect_effect=indirect_effect[-match(intersect(rm_name,rownames(indirect_effect)),rownames(indirect_effect)),]



cairo_pdf("./family/manuscript/draft4/plots/extended/Extended_Figure6.pdf",width=180/25.4, height=95/25.4)
layout(matrix(1:12, 2, 6, byrow=TRUE), widths=c(1.3, 1, 1, 1, 1, 1))
qqplot_adjusted_lambda0(direct_effect$p_all,N_case = 15876,N_control = 15876,title="All (#Trios = 15,876); DE",col_point = "royal blue")
qqplot_adjusted_lambda(direct_effect$p_EUR,N_case = 12813,N_control = 12813,title="EUR (#Trios = 12,813); DE",col_point = "royal blue")
qqplot(direct_effect$p_AFR,title="AFR (#Trios = 792); DE",col_point = "royal blue")
qqplot(direct_effect$p_AMR,title="AMR (#Trios = 1,302); DE",col_point = "royal blue")
qqplot(direct_effect$p_EAS,title="EAS (#Trios = 415); DE",col_point = "royal blue")
qqplot(direct_effect$p_SAS,title="SAS (#Trios = 554); DE",col_point = "royal blue")
qqplot_adjusted_lambda0(indirect_effect$p_all,N_case = 15876,N_control = 15876,title="All (#Trios = 15,876); \u03B4-IDE",col_point = "dark green")
qqplot_adjusted_lambda(indirect_effect$p_EUR,N_case = 12813,N_control = 12813,title="EUR (#Trios = 12,813); \u03B4-IDE",col_point = "dark green")
qqplot(indirect_effect$p_AFR,title="AFR (#Trios = 792); \u03B4-IDE",col_point = "dark green")
qqplot(indirect_effect$p_AMR,title="AMR (#Trios = 1,302); \u03B4-IDE",col_point = "dark green")
qqplot(indirect_effect$p_EAS,title="EAS (#Trios = 415); \u03B4-IDE",col_point = "dark green")
qqplot(indirect_effect$p_SAS,title="SAS (#Trios = 554); \u03B4-IDE",col_point = "dark green")
dev.off()

