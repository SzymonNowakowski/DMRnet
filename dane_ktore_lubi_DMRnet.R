library(devtools)
load_all()
data(insXtr)
data(insytr)

#dmrnet <- DMRnet(Xtr, ytr, family = "gaussian")
#cv.dmrnet  <- cv.DMRnet(Xtr, ytr, family = "gaussian")
#dmr <- DMR(Xtr, ytr, family = "gaussian")
cv.dmr <- cv.DMR(Xtr, ytr, family = "gaussian")
