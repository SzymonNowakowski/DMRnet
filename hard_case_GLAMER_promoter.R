library(devtools)
load_all()
data(promoX)
data(promoy)

############ This problem case is extracted from promoter dataset: cv+sd.glamer failed on the 10th run of massive-promoter.
############ The cause was that grpreg didn't observe the group constraints for the 15th column and returned 2 non-zero betas and one beta equal to 0
############           **** Group constraint explained: https://pbreheny.github.io/grpreg/articles/web/penalties.html#group-selection
############           ****     "These penalties are sparse at the group level – the coefficients within a group will either all equal zero or none will equal zero."
############ Regularization was added to restore the group constraint (either all betas are zero, or all betas are non-zero in a group)



cat("binomial glamer single 4-level factor column:\n")
b_glamer <- DMRnet(X[15], y, family = "binomial", algorithm="glamer")


cat("binomial dmrnet single 4-level factor column:\n")
b_dmrnet <- DMRnet(X[15], y, family = "binomial")

cat("binomial dmr single 4-level factor column:\n")
b_dmr <- DMR(X[15], y, family = "binomial")

cat("binomial glamer:\n")
b_glamer <- DMRnet(X, y, family = "binomial", algorithm="glamer")  #fixed by PR#18-20

cat("binomial dmrnet:\n")
b_dmrnet <- DMRnet(X, y, family = "binomial")

cat("binomial cvg.dmrnet:\n")
b_cv.dmrnet  <- cv.DMRnet(X, y, family = "binomial")

cat("binomial cvg.glamer:\n")
b_cv.dmrnet  <- cv.DMRnet(X, y, family = "binomial", algorithm="glamer")

cat("binomial dmr would not pass, as p>n:\n")
#b_dmr <- DMR(X, y, family = "binomial")

cat("binomial cv.dmr would not pass, as p>n:\n")
#b_cv.dmr <- cv.DMR(X, y, family = "binomial")

y<-(y=="1")


cat("gaussian glamer single 4-level factor column:\n")
glamer <- DMRnet(X[15], y, family = "gaussian", algorithm="glamer")

cat("gaussian dmrnet single 4-level factor column:\n")
dmrnet <- DMRnet(X[15], y, family = "gaussian")


cat("gaussian dmr single 4-level factor column:\n")
dmr <- DMR(X[15], y, family = "gaussian")

cat("gaussian glamer:\n")
glamer <- DMRnet(X, y, family = "gaussian", algorithm = "glamer")  #fixed by PR#18-20

cat("gaussian dmrnet:\n")
dmrnet <- DMRnet(X, y, family = "gaussian")

cat("gaussian cvg.dmrnet:\n")
cv.dmrnet  <- cv.DMRnet(X, y, family = "gaussian")

cat("gaussian cvg.glamer:\n")
cv.dmrnet  <- cv.DMRnet(X, y, family = "gaussian", algorithm="glamer")

cat("gaussian dmr would not pass, as p>n:\n")
#dmr <- DMR(X, y, family = "gaussian")

cat("gaussian cv.dmr would not pass, as p>n:\n")
#cv.dmr <- cv.DMR(X, y, family = "gaussian")




cat("completed\n")
