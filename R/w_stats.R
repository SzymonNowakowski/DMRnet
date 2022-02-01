w_stats <- function(B, VV, ind1, ind2){
    if (ind1 != ind2){
        n <- ind2 - ind1 + 1
        Tmat <- matrix(0, n, n)
        ind <- t(utils::combn(1:n, 2))
        x <- as.matrix(B[ind[, 1]] - B[ind[, 2]])
        #if (ncol(x) == 1) x <- t(x)
        Tmat[ind] <- (x^2)/(VV[cbind(ind[,1],ind[,1])] + VV[cbind(ind[,2],ind[,2])] - 2*VV[cbind(ind[,1],ind[,2])])
        #relating to the 'out[ out<0 ] <- lam' fix from DMR4glm.R and DMR4glm_help.R, it is probably also possible to fix the negative values by changing the line above
        #to the following:   Tmat[ind] <- (x^2)/(VV[cbind(ind[,1],ind[,1])] + VV[cbind(ind[,2],ind[,2])] - VV[cbind(ind[,1],ind[,2])] - VV[cbind(ind[,2],ind[,1])]) that is dropping the assumption VV is symmetrical
        #but this fix is less straightforward and untested.
        #moreover if VV is NOT symmetrical, there is no certainty over what is the value to use it to perform computations, anyway.

        Tmat <- cbind(0, Tmat)
        t_st <- c(0, B^2/VV[cbind(1:n,1:n)])
        Tmat <- rbind(t_st, Tmat)
    } else {
        Tmat <- matrix(c(0, 0, B^2/VV, 0), 2, 2)
    }
    return(Tmat)
}
