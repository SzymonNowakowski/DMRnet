#' DMRnet-package
#'
#' @name DMRnet-package
#' @docType package
#' @description Model selection algorithms for regression and classification, where the predictors can be continuous or categorical and the number of regressors may exceed the number of observations. The selected model consists of a subset of numerical regressors and partitions of levels of factors.
#' @section DMRnet Functions:
#'
#' Similar in use to \pkg{glmnet}. It consists of the following functions:
#'
#' \code{\link{DMR}} - Model selection algorithm for p<n; produces a path of models.
#'
#' \code{\link{DMRnet}} - Model selection algorithm both for p<n and for p>=n; produces a path of models.
#'
#' \code{\link{print.DMR}}, \code{\link{coef.DMR}}, \code{\link{plot.DMR}}, \code{\link{predict.DMR}} - Functions for inspection of the models on the path.
#'
#' \code{\link{gic.DMR}}, \code{\link{cv.DMR}}, \code{\link{cv.DMRnet}} - Functions for final model selection, resulting with one model from the path.
#'
#' \code{\link{coef.gic.DMR}}, \code{\link{coef.cv.DMR}}, \code{\link{plot.gic.DMR}}, \code{\link{plot.cv.DMR}}, \code{\link{predict.gic.DMR}}, \code{\link{predict.cv.DMR}} - Functions for inspection of the final model.
#'
#' \code{miete}, \code{promoter} - Two data sets used for vignettes, examples, etc.
#'
#' For more information see a friendly "Getting started" vignette:
#' @examples
#' \dontrun{
#' vignette("getting-started", package="DMRnet")
#' }
#'
#' @author Agnieszka Prochenka-Sołtys, Piotr Pokarowski, Szymon Nowakowski
#'
#' Maintainer: Szymon Nowakowski <s.nowakowski@mimuw.edu.pl>
#'
#' @references
#'
#' Aleksandra Maj-Kańska, Piotr Pokarowski and Agnieszka Prochenka, 2015. Delete or merge regressors for linear model selection. Electronic Journal of Statistics 9(2): 1749-1778. \doi{10.1214/15-EJS1050}
#'
#' Piotr Pokarowski and Jan Mielniczuk, 2015. Combined l1 and greedy l0 penalized least squares for linear model selection. Journal of Machine Learning Research 16(29): 961-992. \url{https://www.jmlr.org/papers/volume16/pokarowski15a/pokarowski15a.pdf}
#'
#' Piotr Pokarowski, Wojciech Rejchel, Agnieszka Sołtys, Michał Frej and Jan Mielniczuk, 2022. Improving Lasso for model selection and prediction. Scandinavian Journal of Statistics, 49(2): 831–863. \doi{10.1111/sjos.12546}
#'
NULL
