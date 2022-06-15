#' @title coef.gic.DMR
#'
#' @description Extracts coefficients from a \code{gic.DMR} object (for the model with minimal gic).
#'
#' @param object Fitted \code{gic.DMR} object.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Similar to other \code{coef} methods, this function extracts coefficients from a fitted \code{gic.DMR} object for the model with minimal gic.
#'
#' @return Vector of coefficients.
#'
#' @examples
#' data(miete)
#' y <- miete[,1]
#' X <- miete[,-1]
#' m <- DMR(X, y)
#' g <- gic.DMR(m, c = 2.5)
#' coef(g)
#' @export
coef.gic.DMR <- function(object, ...){
        out <- coef.DMR(object$dmr.fit, df = object$df.min)
        return(out[out!=0])
}
