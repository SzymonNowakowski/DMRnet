library(devtools)
load_all()
data(la_svd_Xtr)
data(la_svd_ytr)

############ These problem cases are extracted from insurance dataset:

############## GLAMER with cv+sd failed on 48th run of massive-insurance tests with a seed derived from "insurance_"


  ############ The cause was that grpreg failed with error in Lapack:
  #########Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'
  #########Calls: cv.DMRnet ... <Anonymous> -> newXG -> orthogonalize -> svd -> La.svd



  ######################################### I want to keep it for further evaluation and incpection


cat("glamer:\n")
glamer <- DMRnet(Xtr, ytr, family = "gaussian", algorithm = "glamer")





  ################# ISOLATED CASE ###################
data(crashes_svd)
svd(crashes_svd)    #<--------------crashes



############## GLAMER returned NA model in 3th run of massive-insurance tests with a seed derived from "insUrance_"
data(Glamer_NaN)
mod <- DMRnet(insurance.train.10percent.x, insurance.train.10percent.y, algorithm="glamer")
if (is.na(sum(coef(mod, df=439)))) {    # has NA ?
  stop("Found NA values in a model returned from GLAMER")
}

cat("completed\n")
