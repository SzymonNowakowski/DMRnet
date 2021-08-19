DMRnet4lm_help <- function(S, X, y, fl, clust.method){
    if (sum(S) == 0) {
       mm <- stats::lm.fit(as.matrix(rep(1,length(y))), y)
       return(list(b = c(1, rep(0, sum(fl-1))), rss = sum(mm$res^2)))
    }
    Xn <- X[, S==1, drop = FALSE]
    mfin <- DMR4lm_help(Xn, y, clust.method = clust.method)
    return(mfin)
}
