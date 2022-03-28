library(devtools)
load_all()
data(la_svd_Xtr)
data(la_svd_ytr)

############ This problem case is extracted from insurance dataset: GLAMER with cv+sd failed on 48th run of massive-insurance tests:


############ The cause was that grpreg failed with error in Lapack:
#########Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'
#########Calls: cv.DMRnet ... <Anonymous> -> newXG -> orthogonalize -> svd -> La.svd



######################################### I want to keep it for further evaluation and incpection


cat("glamer:\n")
glamer <- DMRnet(Xtr, ytr, family = "gaussian", algorithm = "glamer")





################# ISOLATED CASE ###################
data(crashes_svd)
svd(crashes_svd)    #<--------------crashes


cat("completed\n")
