t_stats <- function(M, ind1, ind2, sigma_sq, z){
    if (ind1 != ind2){
        n <- nrow(M)
        Tmat <- matrix(0, n, n)
        ind <- t(utils::combn(1:n, 2))
        x <- as.matrix(M[ind[, 1],] - M[ind[, 2], ])
        if (ncol(x) == 1) x <- t(x)
        Tmat[ind] <- ((x%*%z)^2)/(sigma_sq*(apply(x, 1, function(y) t(y)%*%y)))
        Tmat <- cbind(0, Tmat)
        t_st <- c(0, ((M%*%z)^2)/(sigma_sq*(apply(M, 1, function(y) t(y)%*%y))))
        Tmat <- rbind(t_st, Tmat)
    } else {
        Tmat <- matrix(c(0, 0, ((M%*%z)^2)/(sigma_sq*sum(M^2)),0), 2, 2)   #A. Prochenka, PhD Thesis, Algorithm 1, Calculate squared t-statistics for all elementary constraints defined in (2.2)
    }
    return(Tmat)
}
