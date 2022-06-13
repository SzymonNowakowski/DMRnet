#' @title predict.gic.DMR
#'
#' @description Make predictions from a \code{gic.DMR} object (for the model with minimal GIC).
#'
#' @param object Fitted \code{gic.DMR} object.
#'
#' @param newx Data frame of new values for \code{X} at which predictions are to be made.
#'
#' @param type One of: \code{"link"}, \code{"response"}, \code{"class"}. For \code{family="gaussian"} for all values of type it gives the fitted values. For \code{family="binomial"} and \code{type="link"} gives the linear predictors, for \code{type="response"} it gives the fitted probabilities and for \code{type="class"} it produces  the  class  label  corresponding  to  the  maximum  probability.
#'
#' @param unknown.factor.levels The way of handling factor levels in test data not seen while training a model. One of \code{"error"} (the default - throwing an error) or \code{"NA"} (returning \code{NA} in place of legitimate value for problematic rows).
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Similar to other \code{predict} methods, this function predicts fitted values from a fitted \code{gic.DMR} object for the model with minimal GIC.
#'
#' @return Vector of predictions.
#'
#' @examples
#' data(miete)
#' ytr <- miete[1:1500,1]
#' Xtr <- miete[1:1500,-1]
#' Xte <- miete[1501:2053,-1]
#' m <- DMR(Xtr, ytr)
#' g <- gic.DMR(m, c = 2.5)
#' ypr <- predict(g, newx = Xte)
#' @export
predict.gic.DMR <- function(object, newx, type = "link", unknown.factor.levels="error", ...){
        out <- predict.DMR(object$dmr.fit, newx = as.data.frame(newx), df = object$df.min, type = type, unknown.factor.levels=unknown.factor.levels)
        return(out)
}
