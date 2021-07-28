
data(miete)

# miete_shuffled<-data.frame(rent= miete$rent, bathextra = miete$bathextra, year=miete$year, area=miete$area)
#
Xtr <- miete[1:1500,-1]
ytr <- miete[1:1500,1]
Xte <- miete[1501:2053,-1]
#
m1 <- cv.DMRnet(Xtr, ytr)


library(devtools)
load_all()

#m1<-lm.fit(X.te[,1:2000], y.te)

load("~/R/DMRnet/1000rows_5vars_that_cause_problem_for_glmnet.Rdata")
m1<-DMRnet(Xtruncated[,5:1], ytruncated)


subset <- sample(c(1:25000), size = 10)
ssd <- apply(X.tr[subset,], 2, stats::sd)
m1<-cv.DMRnet(X.tr[subset,which(ssd>0)], y.tr[subset])

ssd <- apply(X.tr[subset,], 2, stats::sd)
which(ssd ==0)

#x<-sapply(1:85, function(i) length(which(SS[,i]>0)))
