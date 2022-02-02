glamer_stats <- function(betas_for_all_levels){
  n <- length(betas_for_all_levels)
  #returns a nxn similarity matrix (upper triangular) with 0-diagonal

  Tmat <- matrix(0, n, n)
  for (i in 1:(n-1))
    for (j in (i+1):n)
      Tmat[i,j] <- (betas_for_all_levels[i] - betas_for_all_levels[j])^2

  return(Tmat)
}
