artificial_clustering <- function(points) {
  #points is at least 2-element vector

  # it will mimic hclust results, but it will put all points into the same cluster with height equal to max(abs(points))
 l_infty <- max(abs(points))

 ret <- list()
 ret$height <- rep(l_infty, length(points)-1)
 ret$order_points <- seq_along(points)
 ret$labels <- points
 ret$method <- "artificial"
 ret$dist.method <- "zero_distance"

 ret$merge <- matrix(0, ncol=2, nrow=length(points)-1)
 ret$merge[1, ] <- c(-1, -2)
 for (i in seq_along(points[-c(1, 2)])) {
   ret$merge[i+1, ] <- c(i, -(i+2))
 }
 class(ret) <- "hclust"

 return(ret)
}
