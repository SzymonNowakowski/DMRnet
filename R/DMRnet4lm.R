DMRnet4lm <- function(X, y, clust.method = "complete", o = 5, nlambda = 20, maxp = ceiling(length(y)/2)){
    n <- nrow(X)
    if(n != length(y)){
              stop("Error: non-conforming data: nrow(X) not equal to length(y)")
    }
    ssd <- apply(X, 2, function(x) length(unique(x)))
    if (ssd[1] == 1 & (class(X[,1]) == "numeric" | class(X[,1]) == "integer")){
       X <- X[,-1, drop = FALSE]
       ssd <- ssd[-1]
    }
    if(ncol(X) == 0){
              stop("Error: X has zero columns")
    }
    if(sum(ssd == 1) > 0){
               stop("Error: X has columns with sd = 0 apart from the intercept")
    }
    nn <- sapply(1:ncol(X), function(i) class(X[,i]))
    if(is.null(colnames(X))) colnames(X) <- paste("x", 1:ncol(X), sep = "")
    names(nn) <- colnames(X)
    nn[nn == "integer"] <- "numeric"
    p.x <- ncol(X)
    if(sum(nn != "numeric" & nn != "factor" ) > 0){
              stop("Error: wrong data type, columns should be one of types: integer, factor, numeric")
    }
    faki <- which(nn == "factor")
    n.factors <- length(faki)
    n.levels <- c()
    if (n.factors > 0){
      X[,faki]<-lapply(1:n.factors, function(i) factor(X[,faki[i]]))   #recalculate factors
      n.levels.listed<-sapply(1:n.factors, function(i) levels(X[,faki[i]]))
      n.levels <- sapply(1:n.factors, function(i) length(n.levels.listed[[i]]))
    } else
      n.levels.listed<-c()
    cont <- which(nn == "numeric")
    n.cont <- length(cont)
    ord <- c()
    if(n.cont > 0 ){
              if(n.factors > 0){
                  for (j in 1:n.cont){
                     ord[j] <- sum(n.levels[1:sum(nn[1:cont[j]] == "factor")] - 1) + j + 1
                  }
              } else{
                               ord <- 2:(n.cont + 1)
              }
    }
    X <- X[, order(nn), drop = FALSE]
    nn <- sort(nn)
    x.full <- stats::model.matrix(y~., data = data.frame(y=y, X, check.names = TRUE))
    p <- ncol(x.full)
    fl <- c(n.levels, rep(2, n.cont))
    x.full <- apply(x.full, 2, function(x) sqrt(n/sum(x^2))*x)
    mL <- grpreg::grpreg(x.full[,-1], y, group=rep(1:p.x, fl-1) , penalty = "grLasso", family ="gaussian", nlambda = nlambda)
    RL <- mL$lambda
    dfy <- apply(mL$beta, 2, function(x) sum(x!=0))
    kt <- 1:length(RL)
    if (length(which(dfy >= n)) > 0){
       RL <- RL[-which(dfy >= n)]
       kt <- kt[-which(dfy >= n)]
       dfy <- dfy[-which(dfy >= n)]
    }
    kk <- which(dfy == 0)
    if(length(kk) > 0){
       RL <- RL[-kk]
       kt <- kt[-kk]
    }
    bb <- as.matrix(abs(mL$beta[, kt]))
    SS <- ifelse(bb > 0, 1, 0)
    ii <- duplicated(t(SS))
    prz <- rep(1:p.x, fl-1)
    fac <- apply(bb[-1,ii == FALSE, drop = FALSE], 2, function(x) tapply(x, factor(prz), function(z) sum(z^2)*sqrt(length(z))))
    if(is.null(dim(fac))){
                          fac <- t(as.matrix(fac))
    }
    B <- apply(fac, 2, function(x) stats::quantile(x[x!=0], seq(0, 1, length = (o + 1))[-(o + 1)]))
    B[is.na(B)] <- 0
    S <- sapply(1:o, function(j){
      out <- sapply(1:sum(ii == FALSE), function(i) ifelse(fac[, i] >= B[j,i], 1, 0))
    })
    SS <- matrix(S, p.x, sum(ii == FALSE)*o)
    SS <- t(unique(t(SS)))
    if (p >= n) SS = SS[,-1]
    mm <- lapply(1:ncol(SS), function(i) DMRnet4lm_help(SS[,i], X, y, fl, clust.method))
    maxl <- max(sapply(1:length(mm), function(i) length(mm[[i]]$rss)))
    rss <- sapply(1:length(mm), function(i) c(rep(Inf, maxl - length(mm[[i]]$rss)), mm[[i]]$rss))
    ind <- apply(rss, 1, which.min)
    maxi <- min(p, maxp)
    if (length(ind) > maxi){
       idx <- (length(ind) - maxi):length(ind)
    } else{
       idx <- 1:length(ind)
    }
    be <- sapply(idx, function(i) {b_matrix<-mm[[ind[i]]]$b; if (is.null(dim(b_matrix))) b_matrix<-matrix(b_matrix); part2beta_help(b = b_matrix[,i - sum(rss[, ind[i]] == Inf)], S = SS[,ind[i]], X = X, y = y, fl=fl)})
    #SzN:fixing the above line in accordance with crashes in adult dataset in GLAF
    rownames(be) <- colnames(x.full)
    if(length(ord) > 0){
                  ind1 <- c(1:p)
                  ind1[-ord] = 1:(p - length(ord))
                  ind1[ord] = (p - length(ord) + 1):p
                  be = be[ind1,]
   }
   fit <- list(beta = be, df = length(idx):1, rss = rss[cbind(idx, ind[idx])], n = n, levels.listed = n.levels.listed, arguments = list(family = "gaussian", clust.method = clust.method, o = o, nlambda = nlambda, maxp = maxp), interc = TRUE)
   class(fit) = "DMR"
   return(fit)
}
