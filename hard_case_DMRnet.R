library(devtools)
load_all()
data(insXtr)
data(insytr)

print("dmrnet single 2-level factor column:")
dmrnet <- DMRnet(Xtr[1], ytr, family = "gaussian")
print("dmr single 2-level factor column:")
dmr <- DMR(Xtr[1], ytr, family = "gaussian")

print("dmrnet:")
dmrnet <- DMRnet(Xtr, ytr, family = "gaussian")
print("cv.dmrnet:")
cv.dmrnet  <- cv.DMRnet(Xtr, ytr, family = "gaussian")
print("dmr:")
dmr <- DMR(Xtr, ytr, family = "gaussian")
print("cv.dmr:")
cv.dmr <- cv.DMR(Xtr, ytr, family = "gaussian")

