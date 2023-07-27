library(devtools)
load_all()


############ These problem cases are extracted from insurance dataset:

############ Starting with DMRnet package version 0.3.4.9002 PDMR has functionality as GLAMER in previous versions

############## GLAMER returned NA model in 3rd run of massive-insurance tests with a seed derived from "insUrance_"
cat("PDMR, family gaussian:\n")
load("data/NA_insurance/Glamer_NaN.RData")
mod <- DMRnet(insurance.train.10percent.x, insurance.train.10percent.y, algorithm="PDMR")
if (is.na(sum(coef(mod, df=439)))) {    # has NA ?
  stop("Found NA values in a model returned from GLAMER/gaussian")
}    #FIXED by PR#27

y <- factor(insurance.train.10percent.y>mean(insurance.train.10percent.y))
cat("PDMR, family binomial:\n")
mod <- DMRnet(insurance.train.10percent.x, y, algorithm="PDMR", family="binomial")
if (is.na(sum(coef(mod, df=439)))) {    # has NA ?
  stop("Found NA values in a model returned from GLAMER/binomial")
}


################DMRnet in CVG had problems with NAs too, on the 1st run AFTER the fix in massive_insurance (commit 4eca99d)
#cvg.DMRnet with cvg
#1 error =  NaN -> NaN
#1 df =   ->
#  Error in results_dfs[total_run, index - 1] <- df :
#  replacement has length zero
#Execution halted
#isolated case for DMRnet (CV not necessary here)
cat("DMRnet, family gaussian:\n")
load("data/NA_insurance/DMRnet_NaN.RData")
mod <- DMRnet(X, y)
if (is.na(sum(coef(mod, df=460)))) {    # has NA ?
  stop("Found NA values in a model returned from DMRnet/gaussian")
}     #FIXED by PR#27



cat("DMR, family gaussian:\n")
mod <- DMR(X, y)
if (is.na(sum(coef(mod, df=460)))) {    # has NA ?
  stop("Found NA values in a model returned from DMR/gaussian")
}    #FIXED by PR#28


y <- factor(y>mean(y))
cat("DMRnet, family binomial:\n")
mod <- DMRnet(X, y, family="binomial")
if (is.na(sum(coef(mod, df=460)))) {    # has NA ?
  stop("Found NA values in a model returned from DMRnet/binomial")
}


cat("DMR, family binomial:\n")
mod <- DMR(X, y, family="binomial")
if (is.na(sum(coef(mod, df=460)))) {    # has NA ?
  stop("Found NA values in a model returned from DMR/binomial")
}   #FIXED by PR#28


cat("completed\n")
