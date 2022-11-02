#' @title predict.DMR
#'
#' @description Makes predictions from a \code{DMR} object.
#'
#' @param object Fitted \code{DMR} object.
#'
#' @param newx Data frame of new values for \code{X} at which predictions are to be made. The intercept column should NOT be passed in a call to \code{predict}.
#'
#' @param df Number of parameters in the model for which predictions are required. Default is the entire sequence of models for df=1 to df=p.
#'
#' @param type One of: \code{"link"}, \code{"response"}, \code{"class"}. For \code{family="gaussian"} for all values of \code{type} it gives the fitted values. For \code{family="binomial"} and \code{type="link"} it returns the linear predictors, for \code{type="response"} it returns the fitted probabilities and for \code{type="class"} it produces the class labels corresponding to the maximum probability.
#'
#' @param unknown.factor.levels The way of handling factor levels in test data not seen while training a model. One of \code{"error"} (the default - throwing an error) or \code{"NA"} (returning \code{NA} in place of legitimate value for problematic rows).
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Similar to other \code{predict} methods, this function predicts fitted values from a fitted \code{DMR} object.
#'
#' @return Vector or matrix of predictions.
#'
#' @examples
#' data(miete)
#' ytr <- miete[1:1500,1]
#' Xtr <- miete[1:1500,-1]
#' Xte <- miete[1501:2053,-1]
#' m <- DMR(Xtr, ytr)
#' ypr <- predict(m, newx = Xte, df = 11)
#' @export
predict.DMR <- function(object, newx, df = NULL, type = "link", unknown.factor.levels="error", ...){
         if(is.null(ncol(newx))) stop("Error: newx should be a data frame")

  #  attributing to each column of factors in newx the correct factor list from the original model
         nn <- sapply(1:ncol(newx), function(i) class(newx[,i]))
         factor_columns <- which(nn == "factor")
         n.factors <- length(factor_columns)
         if (n.factors > 0)   #DMRnet only
           if (unknown.factor.levels == "error") {
             for (i in 1:n.factors) {
               newx[,factor_columns[i]] <- factor(newx[,factor_columns[i]])   #start by recalculating the test set factors to minimal possible set
               predict_levels <- levels(newx[,factor_columns[i]])
               if (!min(predict_levels %in% object$levels.listed[[i]])) {#if any factor is outside of the listed levels
                 stop('Error: newx cointains factor levels not known in model training. Replace the unknown factor levels with the known ones and try again or try calling with unknown.factor.levels set to "NA".')
               }
               newx[,factor_columns[i]]<-factor(newx[,factor_columns[i]], levels=object$levels.listed[[i]])   #recalculate factors again, but attribute the original factor list from the train set

             }
           } else if (unknown.factor.levels == "NA") {
             ####first, identify all rows that cause problems in relation to any of the factors
             problematic_rows <- rep(FALSE, nrow(newx))  #start off with all rows set to non-problematic
             for (i in 1:n.factors) {
               train.levels <- object$levels.listed[[i]]
               problematic_rows <- problematic_rows | !(newx[,factor_columns[i]] %in% train.levels)   #factor by factor, update the array of problematic rows
             }
             ####then, remove the problematic rows entirely
             newx<-newx[which(!problematic_rows), , drop=FALSE]
             #### now recalculate factors again, but attribute the original factor list from the train set
             for (i in 1:n.factors) {
               newx[,factor_columns[i]]<-factor(newx[,factor_columns[i]], levels=object$levels.listed[[i]])
             }
           } else stop("Error: wrong unknown.factor.levels: should be one of error, NA")

         dd <- data.frame(newx)
         if (object$interc == TRUE){
            Z <- stats::model.matrix( ~ ., data = dd)
         } else{
            Z <- stats::model.matrix( ~ .-1, data = dd)
         }
         if(ncol(Z) != nrow(object$beta) | is.null(ncol(newx))){
                    stop(paste("Error: non-conforming arrays, newx should be a data frame with ncol equal to", nrow(object$beta)))
         }
         if(!inherits(df, "numeric") & !inherits(df, "integer") & is.null(class(df))){
                      stop("Error: wrong input type, df should have type numeric, integer or NULL")
         }
         if (names(object)[3] == "rss") type = "link"
         if(is.null(df)){
                         out <- Z%*%object$beta
                         colnames(out) <- paste("df", object$df, sep = "")
         } else {
                  if (df>ncol(object$beta)) {
                    stop(paste("Error: requested prediction for model size df =",df,"exceeding maximum model size available which is", ncol(object$beta)))
                  }
                  out <- Z%*%object$beta[,ncol(object$beta) - df + 1]
         }
         if(type == "response" | type == "class"){
                 out <- exp(out)
                 out <- out/(out + 1)
                 if(type == "class"){
                         out <- ifelse(out > .5, 1, 0)
                 }
         } else{
                if (type != "link") stop("Error: wrong type: should be one of link, response, class")
         }


         #### finally, add NA values in place of removed rows
         if (n.factors>0 & unknown.factor.levels == "NA") {
           if (max(problematic_rows)) {    # there actually was at least one problematic_row
             res_out<-rep(NA, length(problematic_rows))
             out_index <- 1
             for (row in 1:length(problematic_rows)) {
               if (!problematic_rows[row]) {
                 res_out[row] <-out[out_index]
                 out_index<-out_index+1
               }
             }
             out<-res_out
           }
         }
         return(out)
}
