glamer_4lm_help <- function(S, betas_with_intercept, mL, X, y, fl, cut_points, clust.method){
  if (sum(S) == 0) {
    mm <- stats::lm.fit(as.matrix(rep(1,length(y))), y)
    return(list(b = c(1, rep(0, sum(fl-1))), rss = sum(mm$res^2), heights = 0, cut_point_model_reference = rep(1, length(cut_points))))
  }

  mfin <- clusters_4lm_help(S, betas_with_intercept, X, y, cut_points, clust.method)
  return(mfin)
}
