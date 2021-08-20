
source("glaf_stats.R")


part2beta_glm_help <- function(b, S, fl){
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
  } else{
    b1 <- b
  }
  return(b1)
}



clusters_4glm_help <- function(S, betas_with_intercept, X, y, clust.method = 'complete', lam = 10^(-7)){

  X <- X[, S==1, drop = FALSE]
  betas_with_intercept <- betas_with_intercept[betas_with_intercept>0]
  betas <- betas_with_intercept[-1]

  n <- nrow(X)
  nn <- sapply(1:ncol(X), function(i) class(X[,i]))
  names(nn) <- colnames(X)
  nn[nn == "integer"] <- "numeric"
  x.full <- stats::model.matrix(y~., data = data.frame(y=y, X, check.names = TRUE))
  p <- ncol(x.full)
  # lmin <- lam*length(y)*2
  # lmax <- lmin*1000
  # RL <- exp(seq(log(lmax), log(lmin), length.out = 20))
  # m <- glmnet::glmnet(x.full, y, lambda = RL, alpha = 0, family = "binomial")
  # be <- c(m$a0[20], m$beta[-1,20])
  faki <- which(nn == "factor")
  n.factors <- length(faki)
  if (n.factors > 0){
    n.levels <- sapply(1:n.factors, function(i) length(levels(X[,faki[i]])))
    p.fac <- sum(n.levels - 1)
  } else{
    p.fac <- 0
  }
  cont <- which(nn == "numeric")
  n.cont <- length(cont)
  namCont <- names(nn)[cont]
  # zb <- exp(x.full%*%be)
  # pix <- zb/(zb + 1)
  # w <- as.numeric(pix*(1-pix))
  # Kan <- t(x.full*w)%*%x.full
  #
  # S <- solve(Kan + diag(rep(2*lam, p)))
  # Var <- S%*%Kan%*%S
  if (n.factors > 0){
    Wmats <- lapply(1:n.factors, function(i) {
      i1 <- ifelse(i == 1, 1, sum(n.levels[1:(i - 1)]-1) +1)
      i2 <- sum(n.levels[1:i]-1)
      out <- glaf_stats(c(0,betas[i1:i2]))   #appending 0 as a beta for the constrained level. Each factor has one level constrained to have beta==0 in AP's PhD Thesis (2.2) and (2.3)
      rownames(out) <- colnames(out) <- levels(X[,faki[i]])
      return(out)
    })
    #cutting dendrograms
    models <- lapply(Wmats, function(x) stats::hclust(stats::as.dist(t(x)), method = clust.method, members = NULL))
    heig <- lapply(1:n.factors, function(x){
      out <- models[[x]]$height
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
  if ((p.fac + 1) < p){  # heights for continous vars are just the betas
    heig.add <- betas[(p.fac + 2):p]^2
    names(heig.add) <- colnames(x.full)[(p.fac + 2):p]
    heig <- c(heig, heig.add)
  }
  heig <- sort(heig)
  len <- length(heig)
  #fitting models on the path
  #Z1 <- Z2 <- c()
  sp <- list()
  form <- c()
  nl <- 0
  Z1 <- X[, faki, drop = FALSE]
  if (n.factors > 0){
    for (i in 1:n.factors){
      sp[[i]] <- 1:n.levels[i]
      sp[[i]][sp[[i]] != 1] <- sp[[i]][sp[[i]] != 1] + nl
      nl <- nl + length(unique(sp[[i]])) - 1
    }
  }
  Z2 <- X[,namCont, drop = FALSE]
  Z <- cbind(Z1,Z2)
  dane <- data.frame(y=y, Z, check.names = T)
  ZZ <- stats::model.matrix(y~., data = dane)

  lmin <- lam*length(y)*2
  lmax <- lmin*1000
  RL <- exp(seq(log(lmax), log(lmin), length.out = 20))

  m <- glmnet::glmnet(ZZ, y, lambda = RL, alpha = 0, family = "binomial")
  b <- c(m$a0[20], m$beta[-1,20])
  names(b) <- colnames(ZZ)
  zb = exp(ZZ%*%betas_with_intercept)  #originally be!!!!!!!!!!
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
        dod <- min(sp[[kt]][sp[[kt]] != 1])
        sp[[kt]] <- stats::cutree(models[[kt]], h = heig[i])
        if(length(sp[[kt]][sp[[kt]] != 1]) > 0){
          sp[[kt]][sp[[kt]] != 1] <- sp[[kt]][sp[[kt]] != 1] + dod - min(sp[[kt]][sp[[kt]] != 1])
        }
        Z1[,kt] <- X[, faki[kt]]
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
      loglik = c(loglik, sum(log(pix)[y == 1]) + sum(log(1-pix)[y == 0]) - lam*sum(m$beta[-1,20]^2))
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
  b = cbind(b, c(m$coef, rep(0, length(heig) - 1)))
  zb = exp(m$coef*rep(1, length(y)))
  pix = zb/(zb + 1)
  loglik = c(loglik, sum(log(pix)[y == 1]) + sum(log(1-pix)[y == 0]))
  return(list(beta = b, loglik = loglik))
}


glaf_4glm_help <- function(S, betas_with_intercept, X, y, fl, clust.method, lam){
  if (sum(S) == 0) {
    m <- stats::glm.fit(as.matrix(rep(1, length(y))), y, family = stats::binomial())
    zb = exp(m$coef*rep(1, length(y)))
    pix = zb/(zb + 1)
    loglik = sum(log(pix)[y == 1]) + sum(log(1-pix)[y == 0])
    b = c(m$coef, rep(0,sum(fl-1)))   #!!!!!!!!!!!important stability change. Added rep(0,...)
    return(list(b = b, loglik = loglik))
  }
  mfin <- clusters_4glm_help(S, betas_with_intercept, X, y, clust.method = clust.method, lam = lam)
  return(mfin)
}



glaf_4glm <- function(X, y, clust.method = "complete", nlambda = 100, lam = 10^(-7), maxp = ceiling(length(y)/4)){
    if (class(y) != "factor"){
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
    if (ssd[1] == 1 & (class(X[,1]) == "numeric" | class(X[,1]) == "integer")){   #removing the intercept from X, as it will be added later to X.full
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
    mL <- grpreg::grpreg(x.full[,-1], y, group=rep(1:p.x, fl-1) , penalty = "grLasso", family ="binomial", nlambda = nlambda)
    RL <- mL$lambda
    dfy <- apply(mL$beta, 2, function(x) sum(x!=0))
    kt <- 1:length(RL)
    lambdas_with_nonzero_beta_number_too_large <- which(dfy >= n)  #(1) removing predictor sets with more predictors than matrix rows
    if (length(lambdas_with_nonzero_beta_number_too_large) > 0){
       RL <- RL[-lambdas_with_nonzero_beta_number_too_large]  #removing them from lambdas
       kt <- kt[-lambdas_with_nonzero_beta_number_too_large]  #and from lambda indices
       dfy <- dfy[-lambdas_with_nonzero_beta_number_too_large]
    }
    lambdas_with_no_betas <- which(dfy == 0)     #(2) removing predictor sets with 0 predictors
    if(length(lambdas_with_no_betas) > 0){
       RL <- RL[-lambdas_with_no_betas]  #removing them from lambdas
       kt <- kt[-lambdas_with_no_betas]  #and from lambda indices
    }
    bb <- as.matrix(abs(mL$beta[, kt]))   #bb is a matrix listing beta values (rows) respective to the net of lambda values (cols)
    bb_predictor_sets <- ifelse(bb > 0, 1, 0)          #bb_predictor_sets is a matrix listing predictor sets (0 or 1 for each predictor)  (rows) respective to the net of lambda values (cols)
    ii <- duplicated(t(bb_predictor_sets))    #detecting duplicated predictor sets
    prz <- rep(1:p.x, fl-1)
    fac <- apply(bb[-1,ii == FALSE, drop = FALSE], 2, function(x) tapply(x, factor(prz), function(z) sum(z^2)*sqrt(length(z))))
      #(3) removing duplicated predictor sets
      #fac is a matrix listing a normalized sum of betas from bb for each predictor values (rows) respective to the net of lambda values (cols)) but without the intercept
            #and without lambdas for which we have duplicated sets of beta
      #fac is representing greek nu from AP PhD Thesis, page 38
      #nrow = #predictors
      #ncol = #active lambdas
#    if(is.null(dim(fac))){      ###this block of code is a mistery to me - removing it for now. Why the transpose ???? If so, maybe a double transpose were the AP's intention
#                         fac <- t(as.matrix(fac))
#    }

        #first, note that sum(ii==FALSE) is a number of predictor sets and it may be smaller than nlambda because of (1), (2), (3)
    SS <- sapply(1:sum(ii == FALSE), function(i) ifelse(fac[, i] > 0, 1, 0))
        #nrow = #predictors
        #ncol = active lambdas
    #if (p >= n) SS <- SS[,-1]    #what is the purpose of that line - it is mistery to me. commenting for now
    bb<-bb[,ii==FALSE, drop=FALSE]   #betas but for active lambdas only
    mm <- lapply(1:ncol(SS), function(i) glaf_4glm_help(SS[,i], bb[,i], X, y, fl, clust.method = clust.method, lam = lam))
    maxl <- max(sapply(1:length(mm), function(i) length(mm[[i]]$loglik)))
    loglik <- sapply(1:length(mm), function(i) c(rep(-Inf, maxl - length(mm[[i]]$loglik)), mm[[i]]$loglik))
    ind <- apply(loglik, 1, which.max)
    maxi <- min(p, maxp)
    if (length(ind) > maxi){
       idx <- (length(ind) - maxi):length(ind)
    } else{
       idx <- 1:length(ind)
    }
    be <- sapply(idx, function(i) {b_matrix<-mm[[ind[i]]]$b; if (is.null(dim(b_matrix))) b_matrix<-matrix(b_matrix); part2beta_glm_help(b = b_matrix[,i - sum(loglik[, ind[i]] == -Inf)], S = SS[,ind[i]], fl=fl)})
    #!!!!!!!!!!!important stability change. Added matrix( ...$b)
    rownames(be) <- colnames(x.full)
    if(length(ord) > 0){
                  ind1 <- c(1:p)
                  ind1[-ord] = 1:(p - length(ord))
                  ind1[ord] = (p - length(ord) + 1):p
                  be = be[ind1,]
   }
   fit <- list(beta = be, df = length(idx):1, loglik = loglik[cbind(idx, ind[idx])], n = n, levels.listed = n.levels.listed ,arguments = list(family = "binomial", clust.method = clust.method, nlambda = nlambda, lam = lam, maxp = maxp), interc = TRUE)
   class(fit) = "DMR"
   return(fit)
}
