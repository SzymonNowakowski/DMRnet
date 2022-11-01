DMR4lm <- function(X, y, clust.method, lam){

    out <- prelasso_common(X, y)
    X <-             out$X                #DMR specific
    factor_columns<- out$factor_columns   #DMR specific
    n <-             out$n
    n.levels <-      out$n.levels   #DMR specific
    n.factors <-     out$n.factors  #DMR specific
    n.cont <-        out$n.cont     #DMR specific
    names.cont <-    out$names.cont #DMR specific
    levels.listed <- out$levels.listed
    #fl <-            out$fl        #not needed in DMR
    x.full <-        out$x.full     #DMR specific
    p <-             out$p
    p.x <-           out$p.x
    p.fac <-         out$p.fac      #DMR specific
    ord <-           out$ord

    #important: x.full (and its successors, like Ro later) is ordered: first intercept, then factors (and their levels), then numeric column

    if (p >= n) {
      stop("Error: p >= n, DMR works only for p < n, use DMRnet instead")
    }
    m <- stats::lm.fit(x.full, y)



      #QR decompostion of the model matrix
    qX <- qr.Q(m$qr, complete=FALSE)   #matrix Q, eq. 3.1 of A.Prochenka PhD Thesis
                                      #explicitly stating that we want partial results (https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/QR.Auxiliaries)
    rX <- qr.R(m$qr) + diag(rep(lam, ncol(x.full))) #matrix R, eq. 3.1 of A.Prochenka PhD Thesis, regularized with a diagonal matrix
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

           spold <- sp[[kt]]
           sp[[kt]] <- stats::cutree(models[[kt]], h = heig[i])
           if(length(sp[[kt]][sp[[kt]] != 1]) > 0){
                                       sp[[kt]][sp[[kt]] != 1] <- sp[[kt]][sp[[kt]] != 1] + min(spold[spold != 1]) - min(sp[[kt]][sp[[kt]] != 1])
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

   min_value <- min(c(abs(m$coef[!is.na(m$coef)]), abs(b[!is.na(b) & (b!=0)])))
   b[is.na(b)] <- min_value / 1000.0
   m$coef[is.na(m$coef)] <- min_value / 1000.0 #setting a very small (close to 0) value for the variables exceeding design matrix rank
                    #consult the comment in part2beta_help() for longer explanation

   b <- cbind(b, c(m$coef, rep(0, length(heig) - 1)))
   rss <- c(rss, sum(m$res^2))

   b <- b[ord,]  #reordering betas to reflect the original matrix X

   fit <- list(beta = b, df = p:1, rss = rss, n = n, levels.listed = levels.listed, lambda=numeric(0), arguments = list(family = "binomial", clust.method = clust.method), interc = TRUE)
   class(fit) = "DMR"
   return(fit)
}
