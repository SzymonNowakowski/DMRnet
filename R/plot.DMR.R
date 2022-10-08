#' @title plot.DMR
#'
#' @description Plots coefficients from a \code{DMR} object.
#'
#' @param x Fitted \code{DMR} object.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Produces a coefficient profile plot of the coefficient paths for a fitted \code{DMR} object.
#'
#' @examples
#' data(miete)
#' y <- miete[,1]
#' X <- miete[,-1]
#' m <- DMR(X, y)
#' plot(m)
#' @export
plot.DMR <- function (x, ...){
         if (x$interc == TRUE){
           graphics::matplot(x$df, t(x$beta[-1,]), type = "l", lty = "solid", ylab = "Coefficients", xlab = "df", ...)
         } else{
           graphics::matplot(x$df, t(x$beta), type = "l", lty = "solid", ylab = "Coefficients", xlab = "df", ...)
         }
}
