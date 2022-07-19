DMRnet4glm <- function(X, y, clust.method, o, nlambda, lam, maxp, lambda) {
    if (!inherits(y, "factor")){
       stop("Error: y should be a factor")
    }
    lev <- levels(factor(y))
    if (length(lev) != 2){
       stop("Error: factor y should have 2 levels")
    }
    y <- ifelse(y == lev[2], 1, 0)
    n <- nrow(X)
    if(n != length(y)){
              stop("Error: non-conforming data: nrow(X) not equal to length(y)")
    }
    ssd <- apply(X, 2, function(x) length(unique(x)))
    if (ssd[1] == 1 & (inherits(X[,1], "numeric") | inherits(X[,1], "integer"))){
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
    factor_columns <- which(nn == "factor")
    n.factors <- length(factor_columns)
    n.levels <- c()
    if (n.factors > 0){
      X[,factor_columns]<-lapply(1:n.factors, function(i) factor(X[,factor_columns[i]]))   #recalculate factors to minimal possible set
      levels.listed<-lapply(1:n.factors, function(i) levels(X[,factor_columns[i]]))
      n.levels <- sapply(1:n.factors, function(i) length(levels.listed[[i]]))
    } else
      levels.listed<-c()
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

    if (is.null(lambda)) {
      user.lambda<-substitute()    #make user.lambda - paradoxically - not present in a call to grpreg
    } else {
      nlambda <- length(lambda)   #override this parameter
      user.lambda <- lambda
    }

    #to keep the x.full as a matrix after the intercept column is removed in case of a single 2-level factor column in X, drop=FALSE must be provided in grpreg call below in x.full[,-1, drop=FALSE]
    #otherwise x.full[,-1] would be reduced to a vector and grpreg would crash
    mL <- grpreg::grpreg(x.full[,-1, drop=FALSE], y, group=rep(1:p.x, fl-1) , penalty = "grLasso", family ="binomial", nlambda = nlambda, lambda = user.lambda)
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
      out <- sapply(1:sum(ii == FALSE), function(i) ifelse(fac[, i] >= B[j,i], 1, 0))#the smallest lambda is the-only-intercept lambda. All other betas are 0. In this way, the o-based quantiles are all-zero and this: 0==fac[, i]>=B[j,i]==0 is TRUE, so a column with value 1 in all rows is returned
      #this is a bit tricky, but this line (and condition `>=` above) indeed determines the first column of SS to be a column with value 1 in all rows (the full model)
    })
    SS <- matrix(S, p.x, sum(ii == FALSE)*o)
    SS <- t(unique(t(SS)))
    if (p >= n) SS = SS[,-1]
    mm <- lapply(1:ncol(SS), function(i) DMRnet4glm_help(SS[,i], X, y, fl, clust.method = clust.method, lam = lam))
    maxl <- max(sapply(1:length(mm), function(i) length(mm[[i]]$loglik)))
    loglik <- sapply(1:length(mm), function(i) c(rep(-Inf, maxl - length(mm[[i]]$loglik)), mm[[i]]$loglik))
    ind <- apply(loglik, 1, which.max)
    maxi <- min(p, maxp)
    if (length(ind) > maxi){
       idx <- (length(ind) - maxi):length(ind)
    } else{
       idx <- 1:length(ind)
    }
    be <- sapply(idx, function(i) {
      #SzN:fixing this in accordance with crashes in adult dataset in GLAMER

      b_matrix<-mm[[ind[i]]]$b;
      if (is.null(dim(b_matrix))) {
        b_matrix<-matrix(b_matrix);  #note this shouldn't be b_matrix<-t(as.matrix(b_matrix)) :   with other matrices that degenerated to HORIZONTAL vectors we want to have HORIZONTAL matrices
        #in this case, however, we want a VERTICAL matrix
        #in b_matrix[,model_index_within_group(i)] we take columns of b_matrix if b_matrix is a legitimate matrix
        #when it is degenerate (a vector), we want that taking the full first column (as model_index_within_group(i) == 1 in those cases) takes this whole vector
      }
      part2beta_glm_help(b = b_matrix[,i - sum(loglik[, ind[i]] == -Inf)], S = SS[,ind[i]], fl=fl)
    })
    rownames(be) <- colnames(x.full)
    if(length(ord) > 0){
                  ind1 <- c(1:p)
                  ind1[-ord] = 1:(p - length(ord))
                  ind1[ord] = (p - length(ord) + 1):p
                  be = be[ind1,]
   }
   fit <- list(beta = be, df = length(idx):1, loglik = loglik[cbind(idx, ind[idx])], n = n, levels.listed = levels.listed, lambda=mL$lambda, arguments = list(family = "binomial", clust.method = clust.method, o = o, nlambda = nlambda, lam = lam, maxp = maxp, lambda = lambda), interc = TRUE)
   class(fit) = "DMR"
   return(fit)
}
