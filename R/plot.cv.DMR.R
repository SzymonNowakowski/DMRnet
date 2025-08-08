#' @title plot.cv.DMR
#'
#' @description Plots cross-validated error values from a \code{cv.DMR} object.
#'
#' @param x Fitted \code{cv.DMR} object.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @details Produces a plot of cross-validated mean error values for the entire sequence of models from the fitted \code{cv.DMR} object. The mean error values are presented with dots. Brackets indicate range of one standard error from the mean error. The \code{df.min} (the smallest model minimizing the mean cross-validation error) and \code{df.1se} (the smallest model falling under minimum mean cross validation error plus its one standard error indicated by the upper bracket) are marked with red and blue dots, respectively.
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
      x_values <- utils::tail(x_values, length(x$cvm))

    point_size <- 1.0

    #establishing the min and max limits of the viewport
    text_height <- graphics::strheight("df.min", cex = 0.7, units="inches") * 1.33 # value of 1.33 includes the shift induced by pos=1)
    viewport_height <- graphics::par("pin")[2]  # height of plot region in inches
    viewport_width <- graphics::par("pin")[1]  # width in inches
    arrowhead_length <- 0.006 * viewport_width  # e.g. 1% of width

    ylim_min <- min(x$cvm - x$cvse)
    ylim_max <- max(x$cvm + x$cvse)
    ylim_min <- ylim_min - (ylim_max - ylim_min) * text_height / viewport_height

    graphics::plot(x_values, x$cvm, pch = 16, cex = point_size, xlab = "df", ylab = "cv error", ylim=c(ylim_min, ylim_max), ...)
    #df.min point & text
    graphics::points(x$df.min, min(stats::na.omit(x$cvm)), pch = 16, col = "red", cex = point_size)
    graphics::text(x$df.min, min(stats::na.omit(x$cvm)), "df.min", pos=1, cex=0.7, col = "red")

    if (!is.null(x$df.1se)) {   # df.1se is defined
      #df.1se point & text
      if (sum(x_values == x$df.1se) == 1) {  # it is present on-screen
        graphics::points(x$df.1se, x$cvm[x_values == x$df.1se], pch = 16, col = "blue", cex = point_size)
        graphics::text(x$df.1se, x$cvm[x_values == x$df.1se], "df.1se", pos=3, cex=0.7, col = "blue")
      }
    }

    # presentation of vertical grey brackets to indicate standard error
    for (i in seq_along(x_values)) {
      se <- x$cvse[i]
      mu <- x$cvm[i]
      if (is.finite(se) && is.finite(mu)) {
        graphics::arrows(x0 = x_values[i], y0 = mu - se,
                         x1 = x_values[i], y1 = mu + se,
                         angle = 90, code = 3, length = arrowhead_length, col = "grey")
      }
    }

  }
