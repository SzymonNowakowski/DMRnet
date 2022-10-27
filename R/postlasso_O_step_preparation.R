postlasso_O_step_preparation <- function(p, p.x, n, o, fac, interc) {

  ncol <- ncol(fac)
  # O - order step, the additional step in SOSnet (hence the name, sOs-net) and DMRnet, but not in glamer

  # called from SOSnet and DMRnet

  B <- apply(fac, 2, function(x) stats::quantile(x[x!=0], seq(0, 1, length = (o + 1))[-(o + 1)]))
  B[is.na(B)] <- 0
  S <- sapply(1:o, function(j){
    out <- sapply(1:ncol, function(i) ifelse(fac[, i] >= B[j,i], 1, 0))
    # the largest (the first) lambda is the-only-intercept lambda. All other betas are 0.
    # even if it failed to be computed this way in glmnet or grpreg, an artificial lambda like this was added in postlasso_common() function

    # In this way, the o-based quantiles are all-zero and this: 0==fac[, i]>=B[j,i]==0 is TRUE, so a column with value 1 in all rows is returned
    # this is a bit tricky, but this line (and condition `>=` above) indeed determines the first column of SS to be a column with value 1 in all rows (the full model)
  })
  SS <- matrix(S, p.x, ncol*o)
  SS <- t(unique(t(SS)))

  if (ncol(SS) == 1) {
    # artificial column(s) must be added - we need at least 2 columns, as one of them may get removed below
    # the first one - the one detected in the if clause above - is the full model,
    if (interc) {   # if the intercept is allowed (cases of DMRnet or SOSnet with interc set to TRUE)
                    # the next one will be the intercept only model
      SS <- cbind(SS, c(rep(0, nrow(SS))))
    } else {   # if the intercept is NOT allowed (case of SOSnet with interc set to FALSE)
               # the next ones will be one-variable models
      SS <- cbind(SS, diag(nrow(SS)))
    }

  }

  if (p >= n) SS = SS[,-1, drop=FALSE]    #this removes the full model in the high-dimensional scenario of p>=n. In the high-dimensional scenario, the full model is too large to be further analysed

  return(SS)
}
