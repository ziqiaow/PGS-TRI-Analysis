#example analysis
set.seed(04262023)

#Generate 200000 families with offsprings disease prevalence to be 1%
dat = sim_prospective_population(n_fam=200000, #Number of families in the population
                               cor_e_prs=F, #Assuming no correlation between PGS and E
                               cor_strat=0.25,
                               rho2=0.25,
                               rho_mf=0, #No assortative mating
                               alpha_fam=-5.5, #Intercept
                               betaG_normPRS=0.4, #Main effect of PGS
                               betaE_bin=0.2, #Main effect of E1
                               betaE_norm=0.2, #Main effect of E2
                               betaGE_normPRS_bin=0, #Interaction effect of PGSxE1
                               betaGE_normPRS_norm=0, #Interaction effect of PGSxE2
                               envir=T)

#Randomly select 1000 affected probands and their families
id = sample(which(dat$D_sim==1),size=1000,replace=F)
PRS_fam=dat$pgs_fam[id,]
envir=dat$E_sim[id,]


startTime <- Sys.time()
res_sim = PGS.TRI(pgs_offspring = PRS_fam[,1], 
                 pgs_mother = PRS_fam[,2], 
                 pgs_father = PRS_fam[,3],
                 GxE_int = TRUE,
                 formula = ~ factor(E_sim_bin)+E_sim_norm,
                 E = envir, 
                 parental_indirect = FALSE,
                 side = 2)
endTime <- Sys.time()
print(endTime - startTime)
#Time difference of 0.005946875 secs

#Print the result
res_sim$res_beta
# Estimate Std.Error   Z.value     Pvalue
# PGS                       0.38144034 0.10852411  3.5147981 0.0004400885
# PGS x factor(E_sim_bin)1 -0.09540741 0.13803817 -0.6911669 0.4894606688
# PGS x E_sim_norm          0.06266047 0.06394153  0.9799652 0.3271033111