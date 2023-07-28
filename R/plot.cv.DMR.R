#' @title plot.cv.DMR
#'
#' @description Plots cross-validated error values from a \code{cv.DMR} object.
#'
#' @param x Fitted \code{cv.DMR} object.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Produces a plot of cross-validated error values for the entire sequence of models from the fitted \code{cv.DMR} object. The horizontal level indicating separation of one standard deviation from the minimum error is indicated with a blue dashed line. The df.min (the smallest model minimizing the cross-validated error) and df.1se (the smallest model falling under the blue dashed line) are marked with red and blue points, respectively.
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
  x_values <- x$dmr.fit$df
  if (length(x_values) > length(x$cvm))
    x_values <- x_values[(length(x_values) - length(x$cvm) + 1):length(x_values)]

  graphics::plot(x_values, x$cvm, pch = 16, xlab = "df", ylab = "cv error", ...)
  #df.min point & text
  graphics::points(x$df.min, min(stats::na.omit(x$cvm)), pch = 16, col = "red")  #min(x$cvm)==x$cvm[length(x$cvm) - x$df.min + 1]
  graphics::text(x$df.min, min(stats::na.omit(x$cvm)), "df.min", pos=3, cex=0.7, col = "red")

  if (!is.null(x$df.1se)) {
    standard_deviation <- stats::sd(stats::na.omit(x$cvm[x$cvm!=Inf & x$cvm!=-Inf]))
    if (standard_deviation != 0) {
      #1 SD line, arrow and text
      y <- min(stats::na.omit(x$cvm)) + standard_deviation
      graphics::lines(c(x_values[1], x_values[length(x_values)]), c(y, y), lty="dashed", col = "blue")
      graphics::arrows(1, min(stats::na.omit(x$cvm)), 1, y, length=0.05, code=3)
      graphics::text(1, y/2+min(stats::na.omit(x$cvm))/2, "1SD", pos=4, cex=0.7)

      #df.1se point & text
      if (sum(x_values == x$df.1se) == 1) {
        graphics::points(x$df.1se, x$cvm[x_values == x$df.1se], pch = 16, col = "blue")
        graphics::text(x$df.1se, x$cvm[x_values == x$df.1se], "df.1se", pos=4, cex=0.7, col = "blue")
      }
    }
  }
}
