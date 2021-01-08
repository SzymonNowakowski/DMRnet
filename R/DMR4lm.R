DMR4lm <- function(X, y, clust.method = 'complete'){
    n <- nrow(X)
    if(is.null(colnames(X))) colnames(X) <- paste("x", 1:ncol(X), sep = "")
    if(n != length(y)){
              stop("Error: non-conforming data: nrow(X) not equal to length(y)")
    }
    ssd <- apply(X, 2, function(x) length(unique(x)))   #number of unique values in each column of X
    if (ssd[1] == 1 & (class(X[,1]) == "numeric" | class(X[,1]) == "integer")){  # removing the first column in case
        # in case is a numeric constant in X
        # i.e. in case it is an Intercept. Other than that, constant columns are NOT allowed
       X <- X[,-1, drop = FALSE]   #drop=FALSE keeps the dimensions of X
       ssd <- ssd[-1]
    }

    if(ncol(X) == 0) {
      stop("Error: X has zero columns")
    }
    if(sum(ssd == 1) > 0) {   #checking if any other constant columns still  exist. if yes -> error
      stop("Error: X has columns with sd = 0 apart from the intercept")
    }
    nn <- sapply(1:ncol(X), function(i) class(X[,i]))
    names(nn) <- colnames(X)
    nn[nn == "integer"] <- "numeric"
    if(sum(nn != "numeric" & nn != "factor" ) > 0) {
      stop("Error: wrong data type, columns should be one of types: integer, factor, numeric")
    }
    x.full <- stats::model.matrix(y~., data = data.frame(y=y, X[,order(nn), drop = FALSE], check.names = TRUE))
      #important: x.full (and its successors, like Ro later) is ordered: first intercept, then factors (and their levels), then numeric column
    p <- ncol(x.full)   # sum over all dimensions of X, e.g. for a factor with k levels it is taking (k-1) as a summand
    if (p >= n) {
       stop("Error: p >= n, DMR works only for p < n, use DMRnet instead")
    }
    m <- stats::lm.fit(x.full, y)
    factor_columns <- which(nn == "factor")
    n.factors <- length(factor_columns)   #number of factors in X
    if (n.factors > 0) {
       n.levels <- sapply(1:n.factors, function(i) length(levels(X[,factor_columns[i]])))
       p.fac <- sum(n.levels - 1)   #sum of factor dimensions. Again, for a factor with k levels it is taking (k-1) as a summand
    } else{
       p.fac <- 0
    }
    cont <- which(nn == "numeric")
    n.cont <- length(cont)   #number of continuous columns in X
    names.cont <- names(nn)[cont]
    ord <- c()
    if (n.cont > 0) {
              if (n.factors > 0) {
                  for (j in 1:n.cont) {
                     ord[j] <- sum(n.levels[1:sum(nn[1:cont[j]] == "factor")] - 1) + j + 1
                  }    #TODO: understand it better. For now, it is a sum of factor dimensions in columns preceding a continuous column, plus a continuous column index plus one (the intercept)
              } else {
                  ord <- 2:(n.cont + 1)   #For each continuous column, it is its index plus one (the intercept)
              }
    }

      #QR decompostion of the model matrix
    qX <- qr.Q(m$qr, complete=FALSE)   #matrix Q, eq. 3.1 of A.Prochenka PhD Thesis
                                      #explicitly stating that we want partial results (https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/QR.Auxiliaries)
    rX <- qr.R(m$qr)   #matrix R, eq. 3.1 of A.Prochenka PhD Thesis
    Ro <- solve(rX)    #matrix R^-1, eq. 3.1 of A.Prochenka PhD Thesis
    #alternatively: z <- qr.qty(m$qr, y)[1:m$qr$rank]   #the rest of the columns are from the complete Q matrix (https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/qr)
    z <- t(qX)%*%y     #vector z, eq. 3.1 of A.Prochenka PhD Thesis
    sigma_sq <- as.numeric((t(m$res)%*%m$res)/(n - p))   #variance estimator sigma hat squared, eq. 3.1 of A.Prochenka PhD Thesis
    #dissimilarity measures - matrices of squared t-statistics for each factor
    if (n.factors > 0){
       Tmats <- lapply(1:n.factors, function(i) {
          i1 <- ifelse(i == 1, 2, sum(n.levels[1:(i - 1)] - 1) + 2)   #first (cumulatively) level of i-th factor
          i2 <- sum(n.levels[1:i] - 1) + 1                            #last (cumulatively) level of i-th factor
                    #in Ro levels come first and numeric columns last
          out <- t_stats(Ro[i1:i2,], ind1 = i1, ind2 = i2, sigma_sq = sigma_sq, z = z)
          rownames(out) <- colnames(out) <- m$xlevels[[i]]
          return(out)
       })
       #cutting dendrograms
       models <- lapply(Tmats, function(x) stats::hclust(stats::as.dist(t(x)), method = clust.method, members = NULL))
       heig <- lapply(1:n.factors, function(x){
            out <- models[[x]]$he
            names(out)<- rep(x, length(out))
            out
       })
       heig <- unlist(heig)
    } else {
        heig <- c()
        models <- list()
    }
    len <- length(heig)
    heig <- c(0,heig)
    names(heig)[1] = "full"
    if ((p.fac + 1) < p){
        if((p.fac + 2) == p){
          heig.add <- ((Ro[(p.fac + 2):p,]%*%z)^2)/(sigma_sq*sum(Ro[(p.fac + 2):p,]^2))
        } else {
          heig.add <- ((Ro[(p.fac + 2):p,]%*%z)^2)/(sigma_sq*(apply(Ro[(p.fac + 2):p, ], 1, function(y) t(y)%*%y)))
        }
        names(heig.add) <- colnames(x.full)[(p.fac + 2):p]
        heig <- c(heig, heig.add)
    }
    heig <- sort(heig)
    len <- length(heig)
    #fitting models on the path
    Z1 <- Z2 <- c()
    sp <- list()
    form <- c()
    nl <- 0
    Z1 <- X[, factor_columns, drop = FALSE]
    if (n.factors > 0){
            for (i in 1:n.factors){
                sp[[i]] <- 1:n.levels[i]
                sp[[i]][sp[[i]] != 1] <- sp[[i]][sp[[i]] != 1] + nl
                nl <- nl + length(unique(sp[[i]])) - 1
            }
     }
     Z2 <- X[,names.cont, drop = FALSE]
     Z <- cbind(Z1,Z2)
     dane <- data.frame(y=y, Z, check.names = T)
     ZZ <- stats::model.matrix(y~., data = dane)
     m <- stats::lm.fit(ZZ, y)
     b <- m$coef
     rss = sum(m$res^2)
     form <- names.cont
     if (len > 2){
       for (i in 2:(len - 1)){
         kt <- names(heig)[i]
         if(length(intersect(kt, names.cont)) > 0){
                          form <- form[-which(form == kt)]
                          Z2 <- Z2[, form, drop = FALSE]

         } else {
           kt <- as.numeric(kt)
           dod <- min(sp[[kt]][sp[[kt]] != 1])
           sp[[kt]] <- stats::cutree(models[[kt]], h = heig[i])
           if(length(sp[[kt]][sp[[kt]] != 1]) > 0){
                                       sp[[kt]][sp[[kt]] != 1] <- sp[[kt]][sp[[kt]] != 1] + dod - min(sp[[kt]][sp[[kt]] != 1])
           }
           Z1[,kt] <- X[, factor_columns[kt]]
           levels(Z1[,kt]) <- sp[[kt]]
           Z1[,kt] <- factor(Z1[,kt])
           if (kt < length(sp)) for( x in (kt+1):length(sp)){ if (length(sp[[x]][sp[[x]]!=1]) > 0 ) sp[[x]][sp[[x]]!= 1] = sp[[x]][sp[[x]]!=1] - 1}
           nl <- nl - 1
         }
         Z <- cbind(Z1[,which(apply(Z1, 2, function(x) length(unique(x))) != 1)], Z2)
         dane <- data.frame(y = y, Z, check.names = T)
         ZZ <- stats::model.matrix(y~., data = dane)
         m <- stats::lm.fit(ZZ, y)
         be <- c(0, m$coef[-1])
         bb <- c()
         if(n.factors > 0){
                       bb <- unlist(sapply(1:length(sp), function(j) sp[[j]][-1]))
         }
         bb2 <- rep(1, n.cont)
         names(bb2) <- names.cont
         if(length(form) > 0){
                    bb2[form] <- (nl + 2):(nl + 1 + length(form))
         }
         bb <- c(bb, bb2)
         b=cbind(b, c(m$coef[1], be[bb]))
         rss = c(rss, sum(m$res^2))
     }
   }
   m <- stats::lm.fit(as.matrix(rep(1, length(y))), y)
   b = cbind(b, c(m$coef, rep(0, length(heig) - 1)))
   rss = c(rss, sum(m$res^2))
   if(length(ord) > 0){
                  ind <- c(1:p)
                  ind[-ord] = 1:(p - length(ord))
                  ind[ord] = (p - length(ord) + 1):p
                  b = b[ind,]
   }
   fit <- list(beta = b, df = p:1, rss = rss, n = n, arguments = list(family = "binomial", clust.method = clust.method), interc = TRUE)
   class(fit) = "DMR"
   return(fit)
}
