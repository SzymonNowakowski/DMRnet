prelasso_binomial <- function(y) {
#called by all 4glm functions:
  # glamer_4glm
  # DMRnet4glm
  # SOSnet4glm
  # DMR4glm

    if (!inherits(y, "factor")){
    stop("Error: y should be a factor")
  }
  lev <- levels(factor(y))
  if (length(lev) != 2){
    stop("Error: factor y should have 2 levels")
  }
  y <- ifelse(y == lev[2], 1, 0)
  return(y)
}
