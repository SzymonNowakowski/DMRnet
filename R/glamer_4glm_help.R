glamer_4glm_help <- function(S, betas_with_intercept, X, y, fl, cut_points, clust.method, lam){
  if (sum(S) == 0) {
    m <- stats::glm.fit(as.matrix(rep(1, length(y))), y, family = stats::binomial())
    zb = exp(m$coef*rep(1, length(y)))
    pix = zb/(zb + 1)
    loglik = sum(log(pix)[y == 1]) + sum(log(1-pix)[y == 0])
    b = c(m$coef, rep(0,sum(fl-1)))   #SzN:that line is in accordance with DMRnet 4lm version and to avoid crashes in adult dataset in GLAMER
    return(list(b = b, loglik = loglik, heights = 0, cut_point_model_reference = rep(1, length(cut_points))))
  }
  mfin <- clusters_4glm_help(S, betas_with_intercept, X, y, cut_points, clust.method, lam)
  return(mfin)
}
