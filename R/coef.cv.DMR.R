#' @title coef.cv.DMR
#'
#' @description Extracts coefficients from a \code{cv.DMR} object (for the model with minimal cross-validated error /the default/ or the smallest model falling under the upper curve of a prediction error plus one standard deviation).
#'
#' @param object Fitted \code{cv.DMR} object.
#'
#' @param md Value of the model dimension parameter at which predictions are required. The default is \code{md="df.min"} value indicating the model minimizing the cross validation error. Alternatively, \code{md="df.1se"} can be used, indicating the smallest model falling under the upper curve of a prediction error plus one standard deviation.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Similar to other \code{coef} methods, this function extracts coefficients from a fitted \code{cv.DMR} object.
#'
#' @return Vector of coefficients.
#'
#' @examples
#' ## cv.DMR for linear regression
#' set.seed(13)
#' data(miete)
#' y <- miete$rent
#' X <- miete$area
#' cv = cv.DMR(X,y)
#' coef(cv)
#'
#' @export
coef.cv.DMR <- function(object, md="df.min", ...){

  if (md=="df.1se" & !is.null(object$df.1se)) {
    out <- coef.DMR(object$dmr.fit, df = object$df.1se)
  } else if (md=="df.1se") {   #object$df.1se is null
    stop("Error: required the smallest model falling under the upper curve of a prediction error plus one standard deviation, but it is not set. Use size=`df.min` instead, for the model minimizing the cross validation error.")
  } else if (md=="df.min") {
    out <- coef.DMR(object$dmr.fit, df = object$df.min)
  } else{
    stop("Error: wrong md, should be one of: df.min, df.1se")
  }

        return(out[out!=0])
}
