SOSnet4glm_help <- function(S, mL, X, y, lam, interc){
    screenPred <- which(S==1)
    s <- sum(S)
    p.x <- ncol(X)
    Z <- as.matrix(X[, screenPred, drop = FALSE])
    if (dim(Z)[2] > 1){
      lmin <- lam*length(y)*2
      lmax <- lmin*1000
      RL <- exp(seq(log(lmax), log(lmin), length.out = 20))
      mm <- glmnet::glmnet(Z, y, lambda = RL, alpha = 0, intercept = interc, family = "binomial")  #SzN per explanation of PP, this is regularized with ridge penalty (alpha=0) to help with computations of singular cases, but not to sparsify the betas as lasso penalty could
      if (interc == FALSE){
         zb <- exp(Z%*%mm$beta[,20])
         pix <- zb/(zb + 1)
         w <- as.numeric(pix*(1-pix))
         #W <- diag(as.numeric(pix*(1-pix)))
         #Kan <- t(Z)%*%W%*%Z
         Kan <- t(Z*w)%*%Z
         S <- solve(Kan + diag(rep(2*lam, ncol(Z))))
         Var <- S%*%Kan%*%S
         T2 <- mm$beta[,20]^2/diag(Var)
      } else{
         Z <- cbind(1, Z)
         bb <- c(mm$a0[20], mm$beta[,20])
         zb <- exp(Z%*%bb)
         pix <- zb/(zb + 1)
         w <- as.numeric(pix*(1-pix))
         #W <- diag(as.numeric(pix*(1-pix)))
         #Kan <- t(Z)%*%W%*%Z
         Kan <- t(Z*w)%*%Z
         S <- solve(Kan + diag(rep(2*lam, ncol(Z))))
         Var <- S%*%Kan%*%S
         T2 <- (bb^2/diag(Var))[-1]
      }
    } else {
       T2 <- 1
    }
    ind <- order(T2, decreasing = T)
    p <- ncol(Z)
    if (interc == FALSE){
      loglikbe <- sapply(1:p, function(i){
           if (i > 1){
           m2 <- glmnet::glmnet(Z[, ind[1:i], drop = FALSE], y, lambda = RL, alpha = 0, intercept = FALSE, family = "binomial")
           zb = exp(Z[,ind[1:i], drop = FALSE]%*%m2$beta[,20])
           pix = zb/(zb + 1)
           loglik = sum(log(pix)[y == 1]) + sum(log(1-pix)[y == 0]) - lam*sum(m2$beta[,20]^2)
           be <- rep(0, p.x)
           be[screenPred[ind[1:i]]] <- m2$beta[,20]
           return(list(loglik = loglik, be = be))
        } else {
           m2 <- stats::glm.fit(Z[,ind[1:i], drop = FALSE], y, family = stats::binomial())
           pix = m2$fitted.values
           loglik = sum(log(pix)[m2$y == 1]) + sum(log(1-pix)[m2$y == 0])
           be <- rep(0, p.x)
           be[screenPred[ind[1:i]]] <- m2$coef
           return(list(loglik = loglik, be = be))
        }
      })
    } else{
      loglikbe <- sapply(0:(p - 1), function(i){
           if (i > 1){
           m2 <- glmnet::glmnet(Z[, ind[1:i] + 1, drop = FALSE], y, lambda = RL, alpha = 0, intercept = TRUE, family = "binomial")
           bb <- c(m2$a0[20], m2$beta[,20])
           zb = exp(cbind(1, Z[, ind[1:i] + 1, drop = FALSE])%*%bb)
           pix = zb/(zb + 1)
           loglik = sum(log(pix)[y == 1]) + sum(log(1-pix)[y == 0]) - lam*sum(bb^2)
           be <- rep(0, p.x + 1)
           be[1] <- bb[1]
           be[screenPred[ind[1:i]] + 1] <- bb[-1]
           return(list(loglik = loglik, be = be))
        } else {
           if(i == 0){
                m2 <- stats::glm.fit(as.matrix(rep(1, length(y))), y, family = stats::binomial())
           } else{
                m2 <- stats::glm.fit(cbind(1, Z[, ind[1] + 1, drop = FALSE]), y, family = stats::binomial())
           }
           pix = m2$fitted.values
           loglik = sum(log(pix)[m2$y == 1]) + sum(log(1-pix)[m2$y == 0])
           be <- rep(0, p.x + 1)
           be[1] <- m2$coef[1]
           if (i == 1){
             be[screenPred[ind[1]] + 1] <- m2$coef[2]
           }
           return(list(loglik = loglik, be = be))
        }
      })
    }
    return(list(loglikbe = loglikbe, ind = screenPred[ind]))
}
