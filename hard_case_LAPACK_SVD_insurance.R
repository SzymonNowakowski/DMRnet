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
cat("isolated_case:\n")
data(crashes_svd)
d<-svd(crashes_svd)    #<--------------crashes



############## GLAMER with cv+sd had problems with SVD routine on the 5th CV run on 17th run of massive-insurance tests with a seed derived from "insUrance_"
cat("glamer:\n")
data("5_cv_MD")
mod <- DMRnet(Xtr, ytr, algorithm="glamer")   #Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'

  ################# ISOLATED CASE ###################
cat("isolated_case:\n")
data(crashes_svd_1)
d<-svd(crashes_svd)    #<--------------crashes



############## GLAMER with cv+sd had problems with SVD routine on the 4th CV run on 7th run of massive-insurance tests with a seed derived from "insURance_" AFTER the fix in massive_insurance (commit 4eca99d)
cat("glamer:\n")
data("4_cv_MD")
mod <- DMRnet(Xtr, ytr, algorithm="glamer")   #Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'

  ################# ISOLATED CASE ###################
cat("isolated_case:\n")
data(crashes_svd_2)
d<-svd(crashes_svd)    #<--------------crashes

cat("completed\n")
