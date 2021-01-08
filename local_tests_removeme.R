
#data(miete)

# miete_shuffled<-data.frame(rent= miete$rent, bathextra = miete$bathextra, year=miete$year, area=miete$area)
#
# Xtr <- miete_shuffled[1:1500,-1]
# ytr <- miete_shuffled[1:1500,1]
# Xte <- miete_shuffled[1501:2053,-1]
#
# m1 <- DMR(Xtr, ytr)


library(devtools)
load_all()

m1<-lm.fit(X.te[,1:2000], y.te)

m1
m1<-DMR(X.tr, y.tr)
