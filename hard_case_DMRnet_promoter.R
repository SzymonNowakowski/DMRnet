library(devtools)
load_all()
data(promoX)
data(promoy)




cat("binomial glamer single 4-level factor column:\n")
b_glamer <- DMRnet(X[15], y, family = "binomial", algorithm="glamer")


cat("binomial dmrnet single 4-level factor column:\n")
b_dmrnet <- DMRnet(X[15], y, family = "binomial")     #fails in 0.2.0 - Error in grpreg: Error: X must be a matrix or able to be coerced to a matrix
                                                         #fixed by PR#5
cat("binomial dmr single 4-level factor column:\n")
b_dmr <- DMR(X[15], y, family = "binomial")           #passes in 0.2.0

cat("binomial glamer:\n")
b_glamer <- DMRnet(X, y, family = "binomial", algorithm="glamer")

cat("binomial dmrnet:\n")
b_dmrnet <- DMRnet(X, y, family = "binomial")        #fails in 0.2.0 - Error in stats::cutree(models[[kt]], h = heig[i]) :
                                                                        #invalid 'tree' ('merge' component)
                                                                        #In addition: Warning messages:
                                                                        #1: In DMR4glm_help(Xn, y, clust.method = clust.method, lam = lam) :
                                                                        #  NAs introduced by coercion
                                                                        #2: In min(sp[[kt]][sp[[kt]] != 1]) :
                                                                        #  no non-missing arguments to min; returning Inf
                                                         #fixed by PR#6

cat("binomial cv.dmrnet:\n")
b_cv.dmrnet  <- cv.DMRnet(X, y, family = "binomial") #fails in 0.2.0 probably because of factor levels conflict between train and test sets - Error: Missing data (NA's) detected in X.  You must eliminate missing data (e.g., by removing cases, removing features, or imputation) before passing X to grpreg
                                                         #fixed by the correct factor level management in CV and
                                                         #fixed by PR#6

cat("binomial dmr:\n")
b_dmr <- DMR(X, y, family = "binomial")              #fails in 0.2.0 - Error in stats::cutree(models[[kt]], h = heig[i]) :
                                                                        #invalid 'tree' ('merge' component)
                                                                        #In addition: Warning messages:
                                                                        #1: In DMR4glm(X, y, clust.method = clust.method, lam = lam) :
                                                                        #  NAs introduced by coercion
                                                                        #2: In min(sp[[kt]][sp[[kt]] != 1]) :
                                                                        #  no non-missing arguments to min; returning Inf

                                                         #fixed by PR#6
cat("binomial cv.dmr:\n")
b_cv.dmr <- cv.DMR(X, y, family = "binomial")        #fails in 0.2.0 probably because of factor levels conflict between train and test sets - Error in stats::hclust(stats::as.dist(t(x)), method = clust.method, members = NULL) :
                                                                                                                                             #  NA/NaN/Inf in foreign function call (arg 10)
                                                         #fixed by the correct factor level management in CV and
                                                         #fixed by PR#6

y<-factor(y=="1")


cat("glamer single 4-level factor column:\n")
glamer <- DMRnet(X[15], y, family = "gaussian", algorithm="glamer")

cat("dmrnet single 4-level factor column:\n")
dmrnet <- DMRnet(X[15], y, family = "gaussian")     #fails in 0.2.0 - Error in grpreg: Error: X must be a matrix or able to be coerced to a matrix
#fixed by PR#5 in 0.3.0
cat("dmr single 4-level factor column:\n")
dmr <- DMR(X[15], y, family = "gaussian")           #passes in 0.2.0

cat("glamer:\n")
glamer <- DMRnet(X, y, family = "gaussian", algorithm = "glamer")

cat("dmrnet:\n")
dmrnet <- DMRnet(X, y, family = "gaussian")        #fails in 0.2.0 - Error in solve.default(rX) :
#system is computationally singular: reciprocal condition number = 8.83141e-33
#fixed by regularizing rX with a very small positive diagonal matrix in DMRnet4lm_help
cat("cv.dmrnet:\n")
cv.dmrnet  <- cv.DMRnet(X, y, family = "gaussian") #fails in 0.2.0 probably because of factor levels conflict between train and test sets - Error: Missing data (NA's) detected in X.  You must eliminate missing data (e.g., by removing cases, removing features, or imputation) before passing X to grpreg
#fixed by the correct factor level management in CV and
#fixed by regularizing rX with a very small positive diagonal matrix in DMRnet4lm_help
cat("dmr:\n")
dmr <- DMR(X, y, family = "gaussian")              #fails in 0.2.0 - Error in solve.default(rX) :
#system is computationally singular: reciprocal condition number = 8.83141e-33
#fixed by regularizing rX with a very small positive diagonal matrix in DMRnet4lm
cat("cv.dmr:\n")
cv.dmr <- cv.DMR(X, y, family = "gaussian")        #fails in 0.2.0 probably because of factor levels conflict between train and test sets - Error in solve.default(rX) :
#Lapack routine dgesv: system is exactly singular: U[433,433] = 0
#fixed by the correct factor level management in CV and
#fixed by regularizing rX with a very small positive diagonal matrix in DMRnet4lm




cat("completed\n")
