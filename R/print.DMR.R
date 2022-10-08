#' @title print.DMR
#'
#' @description Prints a \code{DMR} object.
#'
#' @param x Fitted \code{DMR} object.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Print a summary of the \code{DMR} path at each step along the path.
#'
#' @return The summary is silently returned.
#'
#' @examples
#' data(miete)
#' y <- miete[,1]
#' X <- miete[,-1]
#' m <- DMR(X, y)
#' print(m)
#' @export
print.DMR <- function (x,...){
    cat("\nArguments: ", "\n")
    print(x$arguments)
    cat("Family of models: ", "\n")
    if (names(x)[3] == "rss"){
       print(cbind(Df = x$df, RSS = x$rss))
    } else {
       print(cbind(Df = x$df, loglik = x$loglik))
    }
}
