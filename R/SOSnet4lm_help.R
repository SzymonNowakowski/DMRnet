SOSnet4lm_help <- function(S, mL, X, y, interc){

    screenPred <- which(S == 1)


    if (interc == FALSE){
       QR <- qr(X[, screenPred, drop = FALSE])
    } else{
       X <- cbind(1, X)
       QR <- qr(X[, c(1,screenPred+1), drop = FALSE])
    }
    indices_count <- QR$rank
    W <- qr.R(QR)[1:indices_count, 1:indices_count]  #only rank-by-rank matrix is taken further

    Q<-tryCatch(   #there has been cases of Inf errors: NA/NaN/Inf in foreign function call (arg 1)
              #the reason behind those is numerical instability of QR decomposition when rank is much lower than columns provided
              # a solution is to try this, and if it fails, recompute QR but restricted to pivoted columns within rank
      qr.Q(QR, complete=FALSE)[, 1:indices_count],  #explicitly stating that we want partial results (https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/QR.Auxiliaries)
      error=function(cond) {
        #message("Numerical instability in QR decomposition detected. Will try to handle this. Original error:\n")
        #message(cond)
        #message("\n")
        if (interc == FALSE){  #recompute QR but restricted to pivoted columns within rank
          QR <- qr(X[, screenPred[QR$pivot[1:QR$rank]], drop = FALSE])
        } else{
          QR <- qr(X[, c(1,screenPred+1)[QR$pivot[1:QR$rank]], drop = FALSE])
        }
        return(qr.Q(QR, complete=FALSE)[, 1:indices_count])
      }
    )
    z <- t(Q)%*%y
    bety <- backsolve(W, z)
    W_ <- backsolve(W, diag(indices_count))
    v_bety <- apply(W_^2, 1, sum)
    T2 <- bety^2/v_bety
    if (interc == TRUE) {
      T2 <- T2[-(QR$pivot[1])]   #the index of intercept
      pivot <- QR$pivot[-1] - 1
    } else
      pivot <- QR$pivot
    ind <- sort(T2, decreasing = TRUE, index.return = TRUE)$ix
    #the indices calculated in ind are from pivot permutation. Also, note, that maybe not all indices were taken (only 'QR$rank' indices are taken here)
    if (interc == TRUE){
       QR <- qr(X[, c(1,screenPred[pivot[ind]] + 1), drop = FALSE])
    } else{
       QR <- qr(X[, screenPred[pivot[ind]], drop = FALSE])
    }
    q2 <- (t(qr.Q(QR))%*%y)^2
    rss <- rep(NA, indices_count)
    rss1 <- sum(y^2) - q2[1]
    rss[indices_count] <- rss1
    if (indices_count > 1){
       for (k in (indices_count-1):1){
           rss[k] <- rss[k+1] - q2[indices_count-k+1]
       }
    }
    return(list(rss = rss, ind = screenPred[pivot[ind]]))
}
