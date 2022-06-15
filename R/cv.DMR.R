#' @title cross-validation for DMR
#'
#' @description Does k-fold cross-validation for \code{DMR} and returns a value for df.
#'
#' @param X Input data frame, of dimension n x p; \code{DMR} works only if p<n, for p>=n see \code{\link{DMRnet}}; each row is an observation vector. Columns can be numerical or integer for continuous predictors or factors for categorical predictors.
#'
#' @param y Response variable. Numerical for \code{family="gaussian"} or a factor with two levels for \code{family="binomial"}. For \code{family="binomial"} the last level in alphabetical order is the target class.
#'
#' @param family Response type; one of: \code{"gaussian"}, \code{"binomial"}.
#'
#' @param clust.method Clustering method used for partitioning levels of factors; see function \href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html}{hclust} in package \pkg{stats} for details. \code{clust.method="complete"} is the default.
#'
#' @param lam The amount of penalization in rigde regression (used for logistic regression in order to allow for parameter estimation in linearly separable setups) or the amount of matrix regularization in case of linear regression. Used only for numerical reasons. The default is 1e-7.
#'
#' @param nfolds Number of folds in cross-validation. The default value is 10.
#'
#' @param indexation.mode How the cross validation algorithm should index the models for internal quality comparisons; one of: \code{"GIC"} (the default) for GIC-indexed cross validation, \code{"dimension"}, for model dimension-indexed cross validation.
#'
#' @details \code{cv.DMR} algorithm does cross-validation for \code{DMR} with \code{nfolds} folds. The df for the minimal estimated prediction error is returned.
#'
#' @return An object with S3 class \code{"cv.DMR"} is  returned,  which  is  a  list  with  the  ingredients  of  the  cross-validation fit.
#' \describe{
#'   \item{df.min}{df (number of parameters) of the model with minimal cross-validated error.}
#'   \item{df.1se}{df (number of parameters) of the smallest model falling under the upper curve of a prediction error plus one standard deviation. Only for the indexation.mode equal to \code{"dimension"}, otherwise it is set to \code{NULL}.}
#'   \item{dmr.fit}{Fitted \code{DMR} object for the full data.}
#'   \item{cvm}{The mean cross-validated error for the entire sequence of models.}
#'   \item{foldid}{The fold assignments used.}
#' }
#'
#' @seealso  \code{\link{plot.cv.DMR}} for plotting, \code{\link{coef.cv.DMR}} for extracting coefficients and \code{\link{predict.cv.DMR}} for prediction.
#'
#' @examples
#' ## cv.DMR for linear regression
#' set.seed(13)
#' data(miete)
#' ytr <- miete$rent[1:1500]
#' Xtr <- miete$area[1:1500]
#' Xte <- miete$area[1501:2053]
#' cv <- cv.DMR(Xtr, ytr)
#' print(cv)
#' plot(cv)
#' coef(cv)
#' ypr <- predict(cv, newx = Xte)
#'
#' @export cv.DMR

cv.DMR <- function(X, y, family = "gaussian", clust.method = 'complete', lam = 10^(-7), nfolds = 10, indexation.mode = "GIC"){

        return(cv_indexation.mode_distribute(X, y, nfolds, indexation.mode, DMR, family=family, clust.method=clust.method, lam=lam))   #this way of calling (i.e. var=var) passes the variable names into the ellipsis, otherwise no variable names would be present in the list(...)
}
