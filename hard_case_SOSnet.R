library(digest)
library(devtools)
load_all()


set.seed(strtoi(substr(digest("hard_case_SOSnet", "md5", serialize = FALSE),1,7),16))  #for reproducibility

######################## first put to test identified 5 hard columns ######################################################################
data(thousand_rows_5vars_that_cause_problem_for_glmnet)

n<-nrow(Xtruncated)
Xtruncated_normalized<- apply(Xtruncated, 2, function(x) sqrt(n/sum(x^2))*x)

#############the root cause of problems with this dataset is the inner workings of glmnet, as illustrated by the following examples
#mod<-glmnet::glmnet(Xtruncated, ytruncated, alpha = 1, intercept = TRUE, nlambda = 20, family = "gaussian")   #alpha=1 is pure Lasso
#mod$beta
# the fifth beta is empty - as should be

#mod_normalized<-glmnet::glmnet(Xtruncated_normalized, ytruncated, alpha = 1, intercept = TRUE, nlambda = 20, family = "gaussian")    #alpha=1 is pure Lasso
#mod_normalized$beta
# all are non-zero(SIC !!!!)
#which(Xtruncated_normalized[,4] != Xtruncated_normalized[,5])
# but the fifth and fourth column ARE THE SAME(!)


cat("dmr, 5 cols, family gaussian:\n")
dmr<-DMR(Xtruncated, ytruncated, family = "gaussian")
cat("dmr, 5 cols, family gaussian:\n")
dmr<-DMR(Xtruncated_normalized, ytruncated, family = "gaussian")
cat("dmrnet, 5 cols, family gaussian:\n")
dmrnet<-DMRnet(Xtruncated, ytruncated, nlambda = 20,  family = "gaussian")
cat("dmrnet, 5 cols, family gaussian:\n")
dmrnet<-DMRnet(Xtruncated_normalized, ytruncated, nlambda = 20,  family = "gaussian")



######################## second, put to test the whole dataset ######################################################################
data(data_hard_for_sosnet)



for (i in 1:10) {
  size<-1000
  subset <- sample(c(1:25000), size = size)
  ssd <- apply(X.tr[subset,], 2, stats::sd)
  cat("cv dmrnet, family gaussian, limited rows matrix", size,"x", length(ssd>0), "\n")
  dmrnet<-cv.DMRnet(X.tr[subset,which(ssd>0)], y.tr[subset])

  pred<-predict(dmrnet, X.te[,which(ssd>0)])
}



for (i in 1:10) {
  size<-1000
  subset <- sample(c(1:25000), size = size)
  ssd <- apply(X.tr[subset,], 2, stats::sd)
  colsubset <- sample(c(1:length(which(ssd>0))), size = 0.8 * size)   #10-fold CV makes training sets with roughly 90% of training rows overall, but you can never know exactly
  cat("cv dmr, family gaussian, limited rows and columns matrix", size,"x", 0.8 * size, "\n")
  dmr<-cv.DMR(X.tr[subset,which(ssd>0)[colsubset]], y.tr[subset])

  pred<-predict(dmr, X.te[,which(ssd>0)[colsubset]])
}


y.tr <- factor(y.tr>mean(y.tr))


for (i in 1:10) {
  size<-1000
  subset <- sample(c(1:25000), size = size)
  ssd <- apply(X.tr[subset,], 2, stats::sd)
  cat("cv dmrnet, family binomial, limited rows matrix", size,"x", length(ssd>0), "\n")
  dmrnet<-cv.DMRnet(X.tr[subset,which(ssd>0)], y.tr[subset], family="binomial")
  pred<-predict(dmrnet, X.te[,which(ssd>0)])
}


for (i in 1:10) {
  size<-1000
  subset <- sample(c(1:25000), size = size)
  ssd <- apply(X.tr[subset,], 2, stats::sd)
  colsubset <- sample(c(1:length(which(ssd>0))), size = 0.8 * size)   #10-fold CV makes training sets with roughly 90% of training rows overall, but you can never know exactly
  cat("cv dmr, family binomial, limited rows and columns matrix", size,"x", 0.8 * size, "\n")
  dmr<-cv.DMR(X.tr[subset,which(ssd>0)[colsubset]], y.tr[subset], family="binomial")

  pred<-predict(dmr, X.te[,which(ssd>0)[colsubset]])
}

cat("completed\n")
