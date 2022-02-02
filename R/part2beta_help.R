part2beta_help <- function(b, S, X, y, fl, valid = TRUE) {
  if (!valid) {    #TODO: this valid parameter seems to be related with the cross-points cv, may be removed together with all other related functionality
    return(rep(0, sum(fl-1)+1))
  }

  Z <- data.frame(y)
  if (sum(S == 0) > 0){
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
    b <- b1
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
  be <- c(0, m$coef[-1])
  b[b == 0] <- 1
  be <- be[b]
  be[1] <- m$coef[1]
  return(be)
}
