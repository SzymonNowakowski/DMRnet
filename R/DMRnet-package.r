#' DMRnet-package
#'
#' @name DMRnet-package
#' @docType package
#' @description Model selection algorithms for regression and classification, where the predcitors can be numerical and categorical and the number of regressors exceeds the number of observations. The seleted model consists of a subset of numerical regressors and partitions of levels of factors.
#' @details Similar in use to \pkg{glmnet}. It consists of the following functions:
#'
#' DMR - Model selection algorithm for p<n; produces a path of models.
#'
#' DMRnet - Model slection algorithm even for p>=n; produces a path of models.
#'
#' print.DMR, coef.DMR, plot.DMR, predict.DMR - For inspection of the models on the path.
#'
#' gic.DMR, cv.DMR, cv.DMRnet - For final model selection, resulting with one model from the path.
#'
#' coef.gic.DMR, coef.cv.DMR, plot.gic.DMR, plot.cv.DMR, predict.gic.DMR, predict.cv.DMR - For inspection of the final model.
#'
#' miete, promoter - Two data sets used for methods illustration.
#'
#' @author Agnieszka Prochenka-Sołtys, Piotr Pokarowski
#'
#' Maintainer: Agnieszka Prochenka-Sołtys <ap220756@mimuw.edu.pl>
#'
#' @references
#'
#' Aleksandra Maj-Kańska, Piotr Pokarowski and Agnieszka Prochenka. "Delete or merge regressors for linear model selection." Electronic Journal of Statistics 9.2 (2015): 1749-1778. \url{https://projecteuclid.org/euclid.ejs/1440507392}
#'
#' Piotr Pokarowski and Jan Mielniczuk. "Combined l1 and greedy l0 penalized least squares for linear model selection." Journal of Machine Learning Research 16.5 (2015). \url{http://www.jmlr.org/papers/volume16/pokarowski15a/pokarowski15a.pdf}
#'
#' Agnieszka Prochenka and Piotr Pokarowski. "Delete or Merge Regressors algorithm." 20th European Young Statisticians Meeting, Aug 2017. \url{https://indico.uu.se/event/317/contribution/32/material/paper/0.pdf}
#'
NULL
