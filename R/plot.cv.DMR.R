#' @title plot.cv.DMR
#'
#' @description Plot cross-validated error values from a \code{cv.DMR} object.
#'
#' @param x Fitted \code{cv.DMR} object.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Produces a plot of cross-validated error values for the entire sequence of models from the fitted \code{cv.DMR} object.
#'
#' @examples
#' ## cv.DMR for linear regression
#' set.seed(13)
#' data(miete)
#' y <- miete$rent
#' X <- miete$area
#' cv = cv.DMR(X,y)
#' plot(cv)
#'
#' @export
plot.cv.DMR <- function(x, ...){
  graphics::plot(length(x$cvm):1, x$cvm, pch = 16, xlab = "df", ylab = "cv error", ...)
  graphics::points(x$df.min, min(x$cvm), pch = 16, col = 2)
  #TODO: plot df.1se too, if available (not NULL)
}
