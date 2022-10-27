DMRnet4lm <- function(X, y, clust.method, o, nlambda, lam, maxp, lambda){

    out <- prelasso_common(X, y)
    n <-             out$n
    levels.listed <- out$levels.listed
    fl <-            out$fl
    x.full <-        out$x.full_normalized
    p <-             out$p
    p.x <-           out$p.x
    ord <-           out$ord
    groups <-        out$groups

    if (is.null(lambda)) {
      user.lambda<-substitute()    #make user.lambda - paradoxically - not present in a call to grpreg
    } else {
      nlambda <- length(lambda)   #override this parameter
      user.lambda <- lambda
    }


    mL <-  grpreg::grpreg(x.full[,-1, drop=FALSE], y, group=groups, penalty = "grLasso", family ="gaussian", nlambda = nlambda, lambda = user.lambda)

    bb <-  postlasso_common(mL$lambda, n, mL$beta)
    fac <- postlasso_fac(bb, groups)   #fac must be computed on bb without intercept, it happens internally in postlasso_fac()
    SS <-  postlasso_O_step_preparation(p, p.x, n, o, fac, interc=TRUE)

    mm <-  lapply(1:ncol(SS), function(i) DMRnet4lm_help(SS[,i], X, y, fl, clust.method, lam))  #note that it the original matrix X, not X.full that enters DMRnet4lm_help function

    return(wrap_up_gaussian(mm, p, maxp, SS, fl, X, y, x.full, ord, n, levels.listed, mL, list(family = "gaussian", clust.method = clust.method, o = o, nlambda = nlambda, maxp = maxp, lambda = lambda)))
}
