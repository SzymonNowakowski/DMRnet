glamer_4lm_help <- function(S, betas_with_intercept, X, y, fl, clust.method, lam){
  if (sum(S) == 0) {
    mm <- stats::lm.fit(as.matrix(rep(1,length(y))), y)
    return(list(b = c(1, rep(0, sum(fl-1))), rss = sum(mm$res^2)))
  }

  mfin <- clusters_4lm_help(S, betas_with_intercept, X, y, clust.method, lam)
  return(mfin)
}
