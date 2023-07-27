part2beta_help <- function(b, S, X, y, fl) {

  Z <- data.frame(y)
  if (sum(S == 0) > 0){   #TODO: this whole block of code (with one exception below) seems to be identical to part2beta_glm_help().
                                  #It SHOULD be verified, refactored and renamed (indicating what it really does)
    b1 <- b[1]
    b <- b[-1]
    for (i in 1:length(S)){
      if(S[i] == 1) {
        b1 <- c(b1, b[1:(fl[i] - 1)])
        b <- b[-c(1:(fl[i] - 1))]
      } else{
        b1 <- c(b1, rep(0, (fl[i] - 1)))
      }

    }
    b <- b1     #TODO2: this code is not present in part2beta_glm_help. It may be related to the fact that b is later used again close to the end
  } else{
    b1 <- b
  }

  b1 <- b1[-1]
  for (i in 1:length(fl)){
    if(fl[i] == 2){
      if (b1[1] != 0) {
        Z <- cbind(Z, X[,i])
        colnames(Z)[ncol(Z)] <- colnames(X)[i]
      }
      b1 <- b1[-1]
    } else{
      if(sum(b1[1:(fl[i] - 1)]) != 0){
        Z1 <- X[,i]
        levels(Z1) <- c(0, b1[1:(fl[i] - 1)])
        Z <- cbind(Z, Z1)
        colnames(Z)[ncol(Z)] <- colnames(X)[i]
      }
      b1 <- b1[-c(1:(fl[i] - 1))]

    }
  }
  ZZ <- stats::model.matrix(y~., data = Z)
  m <- stats::lm.fit(ZZ, y)

  m$coef[is.na(m$coef)] <- min(abs(m$coef[!is.na(m$coef)])) / 1000.0   #the problem identified in insurance dataset is that
            #the ZZ matrix may not be full.rank, hence NA result for some variables.
            #as a fix - setting a very small (close to 0) value for the superflous variables
            #this fix takes care of DMRnet calculation for gaussian families, too
            #the corresponding fix was added to DMR4lm and DMR4glm
            #                            and to DMR4glm_help  (for DMRnet/GLAMER for binomial families)

  be <- c(0, m$coef[-1])
  b[b == 0] <- 1
  be <- be[b]
  be[1] <- m$coef[1]

  return(be)
}
