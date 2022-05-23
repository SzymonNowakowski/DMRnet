#' @title cross-validation for DMRnet
#'
#' @description Does k-fold cross-validation for DMR and returns a value for df.
#'
#' @param X Input data frame, of dimension n x p; each row is an observation vector. Columns can be numerical or integer for continuous predictors or factors for categorical predictors.
#'
#' @param y Response variable. Numerical for family="gaussian" or a factor with two levels for family="binomial". For family="binomial" the last level in alphabetical order is the target class.
#'
#' @param family Response type; one of: "gaussian", "binomial".
#'
#' @param clust.method Clustering method used for partitioning levels of factors; see function \href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html}{hclust} in package \pkg{stats} for details.
#'
#' @param o Parameter of the group lasso screening step, described in \code{\link{DMRnet}}.
#'
#' @param nlambda Parameter of the group lasso screening step, described in \code{\link{DMRnet}}.
#'
#' @param lam The amount of penalization in rigde regression (used for logistic regression in order to allow for parameter estimation in linearly separable setups) or the amount of matrix regularization in case of linear regression. Used only for numerical reasons. The default is 1e-7.
#'
#' @param interc Should intercept(s) be fitted (default=TRUE) or set to zero (FALSE). If in X there are any categorical variables, interc=TRUE.
#'
#' @param maxp Maximal number of parameters of the model, smaller values result in quicker computation.
#'
#' @param nfolds Number of folds in cross-validation. The default value is 10.
#'
#' @param indexation.mode How the cross validation algorithm should index the models for internal quality comparisons; one of: "GIC" (the default) for GIC-indexed cross validation, "dimension", for model dimension-indexed cross validation.
#'
#' @param algorithm The algorithm to be used to merge levels; one of: "DMRnet" (the default), "glamer".
#'
#' @details cv.DMRnet algorithm does k-fold cross-validation for DMRnet. The df for the minimal estimated prediction error is returned.
#'
#' @return An object with S3 class "cv.DMR" is  returned,  which  is  a  list  with  the  ingredients  of  the  cross-validation fit.
#' \describe{
#'   \item{df.min}{df (number of parameters) of the model with minimal cross-validated error.}
#'   \item{df.1se}{df (number of parameters) of the smallest model falling under the upper curve of a prediction error plus one standard deviation. Only for the indexation.mode equal to "dimension", otherwise it is set to NULL.}
#'   \item{dmr.fit}{Fitted DMR object for the full data.}
#'   \item{cvm}{The mean cross-validated error for the entire sequence of models.}
#'   \item{foldid}{The fold assignments used.}
#' }
#'
#' @seealso  \code{\link{plot.cv.DMR}} for plotting, \code{\link{coef.cv.DMR}} for extracting coefficients and \code{\link{predict.cv.DMR}} for prediction.
#'
#' @examples
#' ## cv.DMRnet for linear regression
#' set.seed(13)
#' data(miete)
#' ytr <- miete$rent[1:1500]
#' Xtr <- miete$area[1:1500]
#' Xte <- miete$area[1501:2053]
#' cv <- cv.DMRnet(Xtr, ytr)
#' print(cv)
#' plot(cv)
#' coef(cv)
#' ypr <- predict(cv, newx = Xte)
#'
#' @export cv.DMRnet

cv.DMRnet <- function(X, y, family = "gaussian", clust.method = 'complete', o = 5, nlambda = 20, lam = 10^(-7), interc = TRUE, maxp = ifelse(family == "gaussian", ceiling(length(y)/2), ceiling(length(y)/4)), nfolds = 10, indexation.mode = "GIC", algorithm="DMRnet"){

       return(cv_indexation.mode_distribute(X, y, nfolds, indexation.mode, DMRnet, family=family, clust.method=clust.method, o=o, nlambda=nlambda, lam=lam, interc=interc, maxp=maxp, algorithm=algorithm))
        #this way of calling (i.e. var=var) passes the variable names into the ellipsis, otherwise no variable names would be present in the list(...)
}



####################################################################################################################################
#The whole treatment of factors in CV and in DMRnet/predict pair is based on the following analysis:

#The situation is as follows
#Xtr is training data in cross validation or in a regular call via DMRnet->model
#Xte is test data in cross validation or in a regular call via model->predict

#Without loss of generality, let us consider Xtr and Xte to be one column only, with factors.

# A is a true set of all factor levels in Xtr
# B is a true set of all factor levels in Xte
# C=levels(Xtr) is a set of factor levels in original data that Xtr originates from, but it is still assigned to Xtr via the levels() function.
#    As a rule when taking subsets R does not eliminate redundant factors, so C is a superset of A

#There are 4 classes of problems:
# 1. C is a strict superset of A
#    Then if treated naively, DMRnet(...) when constructing a model would throw an error,
#    because we would end up with NaN values in a column dedicated to this superflous factor level (it would happen when a columns gets normalized).
#    The solution to that is very simple. Before the model gets constructed in DMRnet we recalculate the factor level set, C_new. Then C_new=A.
#    SOLVED
# 2. B does not contain a level(s) present in A
#    (sample case: we did sample to Xtr the single Dutch national from the Insurance data set, and he is not present in Xte,
#    because there is only one instance of Dutch national in the whole Insurance data set).
#    As a result predict(...) would throw an error, because expanded model-matrix dimensions would be conflicting.
#    The solution is simple here, too: in constructing a model make a note about true A set (it is stored in levels.listed variable in a model)
#    and then in predict(...) assign the levels of Xte to be equal to A. Only then create the model-matrix.
#    SOLVED
# 3. B contains a factor level(s) not present in A, AND we are doing CV, so we have access to Xtr
#    The solution is to remove the rows with levels that are going to cause problems later in predict(...) from Xte before the prediction
#    The other solution would be to predict using unknown.factor.levels="NA" flag and then eliminate the NAs from comparisons (this one is not used at present)
#    SOLVED
# 4. B contains a factor level(s) not present in A, AND we are NOT doing CV, so we have no access to Xtr.
#    This case is problematic because this situation gets identified too late - we are already in predict(...).
#    At this point, only the model created by DMRnet function and passed to predict(...) is known.
#    We cannot perform inference and we cannot perform any imputation for the problematic data point, either (we don't know Xtr and have no access to it).
#    All that remains is to throw an error (unknown.factor.levels="error", the default) OR
#                            eliminate the problematic rows, predict, and then replenish the result with NAs in place of problematic values (unknown.factor.levels="NA").
#    PROBLEMATIC
