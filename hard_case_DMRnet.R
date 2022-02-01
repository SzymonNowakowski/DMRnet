library(devtools)
load_all()
data(insXtr)
data(insytr)

print("dmrnet single 2-level factor column:")
dmrnet <- DMRnet(Xtr[1], ytr, family = "gaussian")     #fails in 0.2.0 - Error in grpreg: Error: X must be a matrix or able to be coerced to a matrix
                                                       #fixed by PR#5 in 0.3.0
print("dmr single 2-level factor column:")
dmr <- DMR(Xtr[1], ytr, family = "gaussian")           #passes in 0.2.0

print("dmrnet:")
dmrnet <- DMRnet(Xtr, ytr, family = "gaussian")        #fails in 0.2.0 - Error in solve.default(rX) :
                                                                        #system is computationally singular: reciprocal condition number = 8.83141e-33
                                                       #fixed by regularizing rX with a very small positive diagonal matrix in DMRnet4lm_help
print("cv.dmrnet:")
cv.dmrnet  <- cv.DMRnet(Xtr, ytr, family = "gaussian") #fails in 0.2.0 probably because of factor levels conflict between train and test sets - Error: Missing data (NA's) detected in X.  You must eliminate missing data (e.g., by removing cases, removing features, or imputation) before passing X to grpreg
                                                       #fixed by the correct factor level management in CV and
                                                       #fixed by regularizing rX with a very small positive diagonal matrix in DMRnet4lm_help
print("dmr:")
dmr <- DMR(Xtr, ytr, family = "gaussian")              #fails in 0.2.0 - Error in solve.default(rX) :
                                                                        #system is computationally singular: reciprocal condition number = 8.83141e-33
                                                       #fixed by regularizing rX with a very small positive diagonal matrix in DMRnet4lm
print("cv.dmr:")
cv.dmr <- cv.DMR(Xtr, ytr, family = "gaussian")        #fails in 0.2.0 probably because of factor levels conflict between train and test sets - Error in solve.default(rX) :
                                                                                                                                               #Lapack routine dgesv: system is exactly singular: U[433,433] = 0
                                                       #fixed by the correct factor level management in CV and
                                                       #fixed by regularizing rX with a very small positive diagonal matrix in DMRnet4lm



ytr<-factor(ytr>mean(ytr))


print("binomial dmrnet single 2-level factor column:")
b_dmrnet <- DMRnet(Xtr[1], ytr, family = "binomial")     #fails in 0.2.0 - Error in grpreg: Error: X must be a matrix or able to be coerced to a matrix
                                                         #fixed by PR#5
print("binomial dmr single 2-level factor column:")
b_dmr <- DMR(Xtr[1], ytr, family = "binomial")           #passes in 0.2.0

print("binomial dmrnet:")
b_dmrnet <- DMRnet(Xtr, ytr, family = "binomial")        #fails in 0.2.0 - Error in stats::cutree(models[[kt]], h = heig[i]) :
                                                                        #invalid 'tree' ('merge' component)
                                                                        #In addition: Warning messages:
                                                                        #1: In DMR4glm_help(Xn, y, clust.method = clust.method, lam = lam) :
                                                                        #  NAs introduced by coercion
                                                                        #2: In min(sp[[kt]][sp[[kt]] != 1]) :
                                                                        #  no non-missing arguments to min; returning Inf
                                                         #fixed by PR#6

print("binomial cv.dmrnet:")
b_cv.dmrnet  <- cv.DMRnet(Xtr, ytr, family = "binomial") #fails in 0.2.0 probably because of factor levels conflict between train and test sets - Error: Missing data (NA's) detected in X.  You must eliminate missing data (e.g., by removing cases, removing features, or imputation) before passing X to grpreg
                                                         #fixed by the correct factor level management in CV and
                                                         #fixed by PR#6

print("binomial dmr:")
b_dmr <- DMR(Xtr, ytr, family = "binomial")              #fails in 0.2.0 - Error in stats::cutree(models[[kt]], h = heig[i]) :
                                                                        #invalid 'tree' ('merge' component)
                                                                        #In addition: Warning messages:
                                                                        #1: In DMR4glm(X, y, clust.method = clust.method, lam = lam) :
                                                                        #  NAs introduced by coercion
                                                                        #2: In min(sp[[kt]][sp[[kt]] != 1]) :
                                                                        #  no non-missing arguments to min; returning Inf

                                                         #fixed by PR#6
print("binomial cv.dmr:")
b_cv.dmr <- cv.DMR(Xtr, ytr, family = "binomial")        #fails in 0.2.0 probably because of factor levels conflict between train and test sets - Error in stats::hclust(stats::as.dist(t(x)), method = clust.method, members = NULL) :
                                                                                                                                             #  NA/NaN/Inf in foreign function call (arg 10)
                                                         #fixed by the correct factor level management in CV and
                                                         #fixed by PR#6
