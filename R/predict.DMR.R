#' @title predict.DMR
#'
#' @description Make predictions from a DMR object.
#'
#' @param object Fitted DMR object.
#'
#' @param newx Data frame of new values for X at which predictions are to be made.
#'
#' @param df Number of parameters in the model for which predictions are required. Default is the entire sequence of models for df=1 to df=p.
#'
#' @param type One of: link, response, class. For "gaussian" for all values of type it gives the fitted values. For "binomial" type "link" gives the linear predictors, for type "response" it gives the fitted probabilities and for type "class" it produces  the  class  label  corresponding  to  the  maximum  probability.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Similar to other predict methods, this function predicts fitted values from a fitted DMR object.
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
predict.DMR <- function(object, newx, df = NULL, type = "link", ...){
         if(is.null(ncol(newx))){
                                 stop("Error: newx should be a data frame")
         }

  #  attributing to each column of factors in newx the correct factor list from the original model
         nn <- sapply(1:ncol(newx), function(i) class(newx[,i]))
         faki <- which(nn == "factor")
         n.factors <- length(faki)
         if (n.factors > 0)   #DMRnet only
           for (i in 1:n.factors) {
             newx[,faki[i]] <- factor(newx[,faki[i]])   #start by recalculationg the test set factors to minimal possible set
             predict_levels <- levels(newx[,faki[i]])
             if (!min(predict_levels %in% object$levels.listed[[i]])) {#if any factor is outside of the listed levels
               stop("Error: newx cointains factors not known in model calculations. Replace the unknown factors with the known ones and try again.")
             }
             newx[,faki[i]]<-factor(newx[,faki[i]], levels=object$levels.listed[[i]])   #recalculate factors, attribute the original factor list

           }

         dd <- data.frame(newx)
         if (object$interc == TRUE){
            Z <- stats::model.matrix( ~ ., data = dd)
         } else{
            Z <- stats::model.matrix( ~ .-1, data = dd)
         }
         if(ncol(Z) != nrow(object$beta) | is.null(ncol(newx))){
                    stop(paste("Error: non-conforming arrays, newx should be a data frame with ncol equal to", nrow(object$beta)))
         }
         if(class(df) != "numeric" & class(df) != "integer" & is.null(class(df))){
                      stop("Error: wrong input type, df should have type numeric, integer or NULL")
         }
         if (names(object)[3] == "rss") type = "link"
         if(is.null(df)){
                         out <- Z%*%object$beta
                         colnames(out) <- paste("df", object$df, sep = "")
         } else {
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
         return(out)
}
