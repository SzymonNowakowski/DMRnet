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
#' @param indexation.mode How the cross validation algorithm should index the models for internal quality comparisons; one of: "GIC" (the default), "error".
#'
#' @details cv.DMRnet algorithm does k-fold cross-validation for DMRnet. The df for the minimal estimated prediction error is returned.
#'
#' @return An object with S3 class "cv.DMR" is  returned,  which  is  a  list  with  the  ingredients  of  the  cross-validation fit.
#' \describe{
#'   \item{df.min}{df (number of parameters) for the model with minimal cross-validated error.}
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

cv.DMRnet <- function(X, y, family = "gaussian", clust.method = 'complete', o = 5, nlambda = 20, lam = 10^(-7), interc = TRUE, maxp = ifelse(family == "gaussian", ceiling(length(y)/2), ceiling(length(y)/4)), nfolds = 10, indexation.mode = "GIC"){

       return(cv_indexation.mode_distribute(X, y, nfolds, indexation.mode, DMRnet, family, clust.method, o, n_lambda, lam, interc, maxp))
}
