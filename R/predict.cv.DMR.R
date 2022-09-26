#' @title predict.cv.DMR
#'
#' @description Make predictions from a cv.DMR object (for the model with minimal cross-validated error /the default/ or the smallest model falling under the upper curve of a prediction error plus one standard deviation).
#'
#' @param object Fitted cv.DMR object.
#'
#' @param newx Data frame of new values for \code{X} at which predictions are to be made.
#'
#' @param type One of: \code{"link"}, \code{"response"}, \code{"class"}. For \code{family="gaussian"} for all values of \code{type} it gives the fitted values. For \code{family="binomial"} and \code{type="link"} it returns the linear predictors, for \code{type="response"} it returns the fitted probabilities and for \code{type="class"} it produces the class labels corresponding to the maximum probability.
#'
#' @param md Value of the model dimension parameter at which predictions are required. The default is \code{md="df.min"} value indicating the model minimizing the cross validation error. Alternatively, \code{md="df.1se"} can be used, indicating the smallest model falling under the upper curve of a prediction error plus one standard deviation.
#'
#' @param unknown.factor.levels The way of handling factor levels in test data not seen while training a model. One of \code{"error"} (the default - throwing an error) or \code{"NA"} (returning \code{NA} in place of legitimate value for problematic rows).
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Similar to other \code{predict} methods, this function predicts fitted values from a fitted \code{cv.DMR} object.
#'
#' @return Vector of predictions.
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
#' @export
predict.cv.DMR <- function(object, newx, type = "link", md="df.min", unknown.factor.levels="error", ...){
  if (md=="df.1se") {# & !is.null(object$df.1se)) {
    out <- predict.DMR(object$dmr.fit, newx = as.data.frame(newx), df = object$df.1se, type = type, unknown.factor.levels=unknown.factor.levels)
  #} else if (md=="df.1se") {   #object$df.1se is null
  #  stop("Error: required the smallest model falling under the upper curve of a prediction error plus one standard deviation, but it is not set. Use size=`df.min` instead, for the model minimizing the cross validation error.")
  } else if (md=="df.min") {
    out <- predict.DMR(object$dmr.fit, newx = as.data.frame(newx), df = object$df.min, type = type, unknown.factor.levels=unknown.factor.levels)
  } else{
    stop("Error: wrong md, should be one of: df.min, df.1se")
  }

  return(out)
}
