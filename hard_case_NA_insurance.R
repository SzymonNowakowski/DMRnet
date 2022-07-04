library(devtools)
load_all()


############ These problem cases are extracted from insurance dataset:


############## GLAMER returned NA model in 3rd run of massive-insurance tests with a seed derived from "insUrance_"
data(Glamer_NaN)
mod <- DMRnet(insurance.train.10percent.x, insurance.train.10percent.y, algorithm="glamer")
if (is.na(sum(coef(mod, df=439)))) {    # has NA ?
  stop("Found NA values in a model returned from GLAMER")
}


################DMRnet in CVG had problems with NAs too, on the 1st run AFTER the fix in massive_insurance (commit 4eca99d)
#cvg.DMRnet with cvg
#1 error =  NaN -> NaN
#1 df =   ->
#  Error in results_dfs[total_run, index - 1] <- df :
#  replacement has length zero
#Execution halted
#isolated case for DMRnet (CV not necessary here)
data(DMRnet_NaN)
mod <- DMRnet(X, y)
if (is.na(sum(coef(mod, df=460)))) {    # has NA ?
  stop("Found NA values in a model returned from GLAMER")
}


cat("completed\n")
