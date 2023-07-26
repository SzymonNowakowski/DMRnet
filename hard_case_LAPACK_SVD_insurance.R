
############ These problem cases are extracted from insurance dataset:


############## with OpenBLAS which was supposed to be a workaround for LAPACK issue#672
############## I got an error in DMRnet on a matrix that passed SVD without OpenBLAS (cv.DMRnet gaussian test in hard_case_DMRnet_insurance.R test
############## in the 5th CV run
############## the following is the extracted isolated case

cat("DMRnet 5th run CV isolated_case with OpenBLAS:\n")
data(crashes_svd_with_Open_BLAS)
d<-svd(crashes_svd)    #<--------------crashes


library(devtools)
load_all()

############## GLAMER with cv+sd failed on 48th run of massive-insurance tests with a seed derived from "insurance_"


  ############ The cause was that grpreg failed with error in Lapack:
  #########Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'
  #########Calls: cv.DMRnet ... <Anonymous> -> newXG -> orthogonalize -> svd -> La.svd



  ######################################### I want to keep it for further evaluation and incpection


cat("glamer:\n")
load("data/LAPACK_SVD_insurance/la_svd_Xtr.RData")
load("data/LAPACK_SVD_insurance/la_svd_ytr.RData")
glamer <- DMRnet(Xtr, ytr, family = "gaussian", algorithm = "glamer")


  ################# ISOLATED CASE ###################
cat("isolated_case:\n")
data(crashes_svd)
d<-svd(crashes_svd)    #<--------------crashes



############## GLAMER with cv+sd had problems with SVD routine on the 5th CV run on 17th run of massive-insurance tests with a seed derived from "insUrance_"
cat("glamer:\n")
load("data/LAPACK_SVD_insurance/5_cv_MD.RData")
mod <- DMRnet(Xtr, ytr, algorithm="glamer")   #Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'

  ################# ISOLATED CASE ###################
cat("isolated_case:\n")
data(crashes_svd_1)
d<-svd(crashes_svd)    #<--------------crashes



############## GLAMER with cv+sd had problems with SVD routine on the 4th CV run on 7th run of massive-insurance tests with a seed derived from "insURance_" AFTER the fix in massive_insurance (commit 4eca99d)
cat("glamer:\n")
load("data/LAPACK_SVD_insurance/4_cv_MD.RData")
mod <- DMRnet(Xtr, ytr, algorithm="glamer")   #Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'

  ################# ISOLATED CASE ###################
cat("isolated_case:\n")
data(crashes_svd_2)
d<-svd(crashes_svd)    #<--------------crashes

cat("completed\n")
