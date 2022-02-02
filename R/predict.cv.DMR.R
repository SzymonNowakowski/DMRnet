#' @title predict.cv.DMR
#'
#' @description Make predictions from a cv.DMR object (for the model with minimal cross-validated error).
#'
#' @param object Fitted cv.DMR object.
#'
#' @param newx Data frame of new values for X at which predictions are to be made.
#'
#' @param md Value of the model dimension parameter at which predictions are required. The default is md="df.min" value indicating the model minimizing the cross validation error. Alternatively, md="df.1se" can be used, indicating the smallest model falling under the upper curve of a prediction error plus one standard deviation.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @param type One of: link, response, class. For "gaussian" for all values of type it gives the fitted values. For "binomial" type "link" gives the linear predictors, for type "response" it gives the fitted probabilities and for type "class" it produces  the  class  label  corresponding  to  the  maximum  probability.
#'
#' @details Similar to other predict methods, this function predicts fitted values from a fitted cv.DMR object for the model with minimal cross-validated error.
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
predict.cv.DMR <- function(object, newx, type = "link", md="df.min", ...){
  if (md=="df.1se" & !is.null(object$df.1se)) {
    out <- predict.DMR(object$dmr.fit, newx = as.data.frame(newx), df = object$df.1se, type = type)
  } else if (md=="df.1se") {   #object$df.1se is null
    stop("Error: required the smallest model falling under the upper curve of a prediction error plus one standard deviation, but it is not set. Use size=`df.min` instead, for the model minimizing the cross validation error.")
  } else if (md=="df.min") {
    out <- predict.DMR(object$dmr.fit, newx = as.data.frame(newx), df = object$df.min, type = type)
  } else{
    stop("Error: wrong md, should be one of: df.min, df.1se")
  }

  return(out)
}
