DMR4glm <- function(X, y, clust.method, lam){

    y <- prelasso_binomial(y)

    out <- prelasso_common(X, y)
    X <-             out$X                #DMR specific
    factor_columns<- out$factor_columns   #DMR specific
    n <-             out$n
    n.levels <-      out$n.levels   #DMR specific
    n.factors <-     out$n.factors  #DMR specific
    n.cont <-        out$n.cont     #DMR specific
    namCont <-       out$names.cont #DMR specific
    levels.listed <- out$levels.listed
    #fl <-            out$fl        #not needed in DMR
    x.full <-        out$x.full     #DMR specific
    p <-             out$p
    p.x <-           out$p.x
    p.fac <-         out$p.fac      #DMR specific
    ord <-           out$ord

    #important: x.full (and its successors, like Ro later) is ordered: first intercept, then factors (and their levels), then numeric column

    if (p >= n){
        stop("Error: p >= n, DMR works only for p < n, use DMRnet instead")
    }

    lmin <- lam*length(y)*2
    lmax <- lmin*1000
    RL <- exp(seq(log(lmax), log(lmin), length.out = 20))
###############LASSO#####################
    m <- glmnet::glmnet(x.full, y, lambda = RL, alpha = 0, family = "binomial")  #SzN per explanation of PP, this is regularized with ridge penalty (alpha=0) to help with computations of singular cases, but not to sparsify the betas as lasso penalty could
###########################################
    be <- c(m$a0[20], m$beta[-1,20])

    zb <- exp(x.full%*%be)
    pix <- zb/(zb + 1)
    w <- as.numeric(pix*(1-pix))
    Kan <- t(x.full*w)%*%x.full
    #W <- diag(as.numeric(pix*(1-pix)))
    #Kan <- t(x.full)%*%W%*%x.full
    S <- solve(Kan + diag(rep(2*lam, p)))
    Var <- S%*%Kan%*%S
    if (n.factors > 0){
       Wmats <- lapply(1:n.factors, function(i) {
          i1 <- ifelse(i == 1, 2, sum(n.levels[1:(i - 1)] - 1) + 2)
          i2 <- sum(n.levels[1:i] - 1) + 1
          out <- w_stats(be[i1:i2], Var[i1:i2, i1:i2], ind1 = i1, ind2 = i2)

          out[ out<0 ] <- lam    #this fix is for DMR, it replaces negative values with a very small positive number. The reason negative values are out there is numerical instability when Kan is very close to 0 and Var is not symmetric, then w_stats produces negative numbers

          rownames(out) <- colnames(out) <- levels(X[,factor_columns[i]])
          return(out)
       })
       #cutting dendrograms
       models <- lapply(Wmats, function(x) stats::hclust(stats::as.dist(t(x)), method = clust.method, members = NULL))
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
        heig.add <- be[(p.fac + 2):p]^2/Var[cbind((p.fac + 2):p,(p.fac + 2):p)]
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
     Z2 <- X[, namCont, drop = FALSE]
     Z <- cbind(Z1,Z2)
     dane <- data.frame(y=y, Z, check.names = T)
     ZZ <- stats::model.matrix(y~., data = dane)
     m <- glmnet::glmnet(ZZ, y, lambda = RL, alpha = 0, family = "binomial")
     b <- c(m$a0[20], m$beta[-1,20])
     names(b) <- colnames(ZZ)
     zb = exp(ZZ%*%be)
     pix = zb/(zb + 1)
     loglik = sum(log(pix)[y == 1]) + sum(log(1-pix)[y == 0]) - lam*sum(m$beta@x^2)
     form <- namCont
     if (len > 2){
       for (i in 2:(len - 1)){
         kt <- names(heig)[i]
         if(length(intersect(kt, namCont)) > 0){
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
         m <- glmnet::glmnet(ZZ, y, lambda = RL, alpha = 0, family = "binomial")
         be <- c(m$a0[20], m$beta[-1,20])
         zb = exp(ZZ%*%be)
         pix = zb/(zb + 1)
         loglik = c(loglik, sum(log(pix)[y == 1]) + sum(log(1-pix)[y == 0]) - lam*sum(m$beta@x^2))
         be[1] <- 0
         bb <- c()
         if(n.factors > 0){
                       bb <- unlist(sapply(1:length(sp), function(j) sp[[j]][-1]))
         }
         bb2 <- rep(1, n.cont)
         names(bb2) <- namCont
         if(length(form) > 0){
                    bb2[form] <- (nl + 2):(nl + 1 + length(form))
         }
         bb <- c(bb, bb2)
         b=cbind(b, c(m$a0[20], be[bb]))
     }
   }
   m <- stats::glm.fit(as.matrix(rep(1, length(y))), y, family = stats::binomial())

   min_value <- min(c(abs(m$coef[!is.na(m$coef)]), abs(b[!is.na(b) & (b!=0)])))
   b[is.na(b)] <- min_value / 1000.0
   m$coef[is.na(m$coef)] <- min_value / 1000.0 #setting a very small (close to 0) value for the variables exceeding design matrix rank
                            #consult the comment in part2beta_help() for longer explanation

   b = cbind(b, c(m$coef, rep(0, length(heig) - 1)))
   zb = exp(m$coef*rep(1, length(y)))
   pix = zb/(zb + 1)
   loglik = c(loglik, sum(log(pix)[y == 1]) + sum(log(1-pix)[y == 0]))
   if(length(ord) > 0){
                  ind <- c(1:p)
                  ind[-ord] = 1:(p - length(ord))
                  ind[ord] = (p - length(ord) + 1):p
                  b = b[ind,]
   }
   fit <- list(beta = b, df = p:1, loglik = loglik, n = n, levels.listed = levels.listed, lambda=numeric(0), arguments = list(family = "binomial", clust.method = clust.method, lam = lam), interc = TRUE)
   class(fit) = "DMR"
   return(fit)
}
