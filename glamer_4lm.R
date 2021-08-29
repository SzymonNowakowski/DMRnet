source("glamer_stats.R")

part2beta_4lm_help <- function(b, S, X, y, fl){
  Z <- data.frame(y)
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
    b <- b1
  } else{
    b1 <- b
  }

  b1 <- b1[-1]
  for (i in 1:length(fl)){
    if(fl[i] == 2){
      if (b1[1] != 0) {
        Z <- cbind(Z, X[,i])
        colnames(Z)[ncol(Z)] <- colnames(X)[i]
      }
      b1 <- b1[-1]
    } else{
      if(sum(b1[1:(fl[i] - 1)]) != 0){
        Z1 <- X[,i]
        levels(Z1) <- c(0, b1[1:(fl[i] - 1)])
        Z <- cbind(Z, Z1)
        colnames(Z)[ncol(Z)] <- colnames(X)[i]
      }
      b1 <- b1[-c(1:(fl[i] - 1))]

    }
  }
  ZZ <- stats::model.matrix(y~., data = Z)
  m <- stats::lm.fit(ZZ, y)
  be <- c(0, m$coef[-1])
  b[b == 0] = 1
  be <- be[b]
  be[1] <- m$coef[1]
  return(be)
}


clusters_4lm_help <- function(S, betas_with_intercept, X, y, cut_points, clust.method){

  X <- X[, S==1, drop = FALSE]
  betas_with_intercept <- betas_with_intercept[betas_with_intercept>0]
  betas <- betas_with_intercept[-1]

  n <- nrow(X)
  nn <- sapply(1:ncol(X), function(i) class(X[,i]))
  names(nn) <- colnames(X)
  nn[nn == "integer"] <- "numeric"
  x.full <- stats::model.matrix(y~., data = data.frame(y=y, X, check.names = TRUE))
  p <- ncol(x.full)
  #m <- stats::lm.fit(x.full, y)
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
  #QR decompostion of the model matrix
  # qX <- qr.Q(m$qr)
  # rX <- qr.R(m$qr)
  # Ro <- solve(rX)
  # z <- t(qX)%*%y
  # sigma <- as.numeric((t(m$res)%*%m$res)/(n - p))
  #dissimilarity measures - matrices of squared t-statistics for each factor
  if (n.factors > 0){
    Tmats <- lapply(1:n.factors, function(i) {
      i1 <- ifelse(i == 1, 1, sum(n.levels[1:(i - 1)]-1) +1)
      i2 <- sum(n.levels[1:i]-1)
      out <- glamer_stats(c(0,betas[i1:i2]))   #appending 0 as a beta for the constrained level. Each factor has one level constrained to have beta==0 in AP's PhD Thesis (2.2) and (2.3)
      rownames(out) <- colnames(out) <- levels(X[,faki[i]])
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
    # if((p.fac + 2) == p){
    #   heig.add <- ((Ro[(p.fac + 2):p,]%*%z)^2)/(sigma*sum(Ro[(p.fac + 2):p,]^2))
    # } else {
    #   heig.add <- ((Ro[(p.fac + 2):p,]%*%z)^2)/(sigma*(apply(Ro[(p.fac + 2):p, ], 1, function(y) t(y)%*%y)))
    # }
    heig.add <- betas_with_intercept[(p.fac + 2):p]^2
    names(heig.add) <- colnames(x.full)[(p.fac + 2):p]
    heig <- c(heig, heig.add)
  }
  heig <- sort(heig)
  len <- length(heig)
  #fitting models on the path
  sp <- list()
  form <- c()
  nl <- 0
  if (n.factors > 0){
    for (i in 1:n.factors){
      sp[[i]] <- 1:n.levels[i]
      sp[[i]][sp[[i]] != 1] <- sp[[i]][sp[[i]] != 1] + nl
      nl <- nl + length(unique(sp[[i]])) - 1
    }
  }
  b <- 1:p
  names(b) <- colnames(x.full)
  A <- c()
  form <- namCont

  cut_point_model_reference <- rep(NA, length(cut_points)) #initiated as pointing to the no-model
        #or numeric(0) if the cut_points==NULL
    ###only 1st cut point should reference to the 1st model
  cut_point_model_reference[1] <- 1
  if (len >= 2){
    for (i in 2:(len)){
      a <- rep(0,p)
      kt <- names(heig)[i]
      if (length(cut_points) > 0) {
        ###all cut points located between this heig[i] and the previous heig[i-1] should back-reference to the model we create in this step
        cut_point_model_reference[cut_points > heig[i-1] & cut_points <= heig[i]] <- i
      }
      if(length(intersect(kt, namCont)) > 0){
        jj <- which(form == kt)
        form <- form[-which(form == kt)]
        jj <- which(namCont == kt)
        a[p.fac + jj + 1] <- 1

      } else {
        kt <- as.numeric(kt)
        dod <- min(sp[[kt]][sp[[kt]] != 1])
        spold <- sp[[kt]]
        sp[[kt]] <- stats::cutree(models[[kt]], h = heig[i])
        if(length(sp[[kt]][sp[[kt]] != 1]) > 0){
          sp[[kt]][sp[[kt]] != 1] <- sp[[kt]][sp[[kt]] != 1] + dod - min(sp[[kt]][sp[[kt]] != 1])
        }
        ii <- min(which(spold != sp[[kt]]))
        suma <- ifelse(kt == 1, 0, sum(n.levels[1:(kt-1)] - 1))
        if(sp[[kt]][ii] == 1){
          a[suma + ii] <- 1
        } else {
          a[suma + ii] <- 1
          a[suma + min(which(sp[[kt]] == sp[[kt]][ii]))] <- -1
        }
        if (kt < length(sp)) for( x in (kt+1):length(sp)){ if (length(sp[[x]][sp[[x]]!=1]) > 0 ) sp[[x]][sp[[x]]!= 1] = sp[[x]][sp[[x]]!=1] - 1}
        nl <- nl - 1
      }
      A <- cbind(A, a)
      be <- c(0, 2:(p-i+2))
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
      b=cbind(b, c(1, be[bb]))
    }
  }
  A[1,] <- rep(0, p-1)
  m <- stats::lm.fit(x.full, y)  #moved here from top
  rX <- qr.R(m$qr)               #moved here from top
  qX <- qr.Q(m$qr)               #moved here from top
  z <- t(qX)%*%y                 #moved here from top
  S <- forwardsolve(t(rX), A)
  QRs <- qr(S)
  W <- qr.Q(QRs)
  wyn <- (t(W)%*%z)^2
  len <- nrow(wyn)
  Tr <- round(lower.tri(matrix(1, len, len))) + diag(rep(1, len))
  r22 <- Tr%*%wyn
  RSS <- (sum(y^2) - sum(z^2))
  RSS2 <- c(RSS, as.vector(RSS + r22))
  return(list(b = b, rss = RSS2, heights = heig, cut_point_model_reference = cut_point_model_reference))
}


glamer_4lm_help <- function(S, betas_with_intercept, mL, X, y, fl, cut_points, clust.method){
  if (sum(S) == 0) {
    mm <- stats::lm.fit(as.matrix(rep(1,length(y))), y)
    return(list(b = c(1, rep(0, sum(fl-1))), rss = sum(mm$res^2), heights = 0, cut_point_model_reference = rep(1, length(cut_points))))
  }

  mfin <- clusters_4lm_help(S, betas_with_intercept, X, y, cut_points, clust.method)
  return(mfin)
}

glamer_4lm <- function(X, y, cut_points = NULL, clust.method = "complete", nlambda = 100, maxp = ceiling(length(y)/2)){
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
      n.levels.listed<-lapply(1:n.factors, function(i) levels(X[,faki[i]]))
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
    bb <- as.matrix(abs(mL$beta[, kt])) #bb is a matrix listing beta values (rows) respective to the net of lambda values (cols)
    bb_predictor_sets <- ifelse(bb > 0, 1, 0)          #bb_predictor_sets is a matrix listing predictor sets (0 or 1 for each predictor)  (rows) respective to the net of lambda values (cols)
    ii <- duplicated(t(bb_predictor_sets))    #detecting duplicated predictor sets
    prz <- rep(1:p.x, fl-1)
    fac <- apply(bb[-1,ii == FALSE, drop = FALSE], 2, function(x) tapply(x, factor(prz), function(z) sum(z^2)*sqrt(length(z))))
    # if(is.null(dim(fac))){
    #                       fac <- t(as.matrix(fac))
    # }

    SS <- sapply(1:sum(ii == FALSE), function(i) ifelse(fac[, i] > 0, 1, 0))

    #if (p >= n) SS = SS[,-1]

    bb<-bb[,ii==FALSE, drop=FALSE]   #betas but for active lambdas only

    mm <- lapply(1:ncol(SS), function(i) glamer_4lm_help(SS[,i], bb[,i], mL, X, y, fl, cut_points, clust.method))
    maxl <- max(sapply(1:length(mm), function(i) length(mm[[i]]$rss)))
    rss <- sapply(1:length(mm), function(i) c(rep(Inf, maxl - length(mm[[i]]$rss)), mm[[i]]$rss))

    if (length(cut_points)>0) {
      cut_point_model_reference <- sapply(1:length(mm), function(i) mm[[i]]$cut_point_model_reference)
      mask <- matrix(Inf, nrow=nrow(rss), ncol=ncol(rss))   #it will mask by Inf all models not back-referenced in cut_points
      for (i in 1:length(mm))
        mask[maxl - na.omit(cut_point_model_reference[,i]) + 1, i] <- 0
      rss_with_mask <- rss + mask
    } else
      rss_with_mask <- rss
    ind <- apply(rss_with_mask, 1, which.min)  #in each row, which is the index of a model minimizing rss
    indInf<- apply(rss_with_mask,1,min) == Inf

    maxi <- min(p, maxp)
    if (length(ind) > maxi){
       idx <- (length(ind) - maxi):length(ind)   #but real model sizes are still length(idx):1
    } else{
       idx <- 1:length(ind)
    }

    #smallest models are last
    model_group <- function(i) {ind[i]}
    model_index_within_group <- function(i) {i - sum(rss[, model_group(i)] == Inf)}
    model_index_within_group_inverted <- function(i) {length(mm[[model_group(i)]]$heights)-model_index_within_group(i)+1}

    be <- sapply(idx, function(i) {b_matrix<-mm[[model_group(i)]]$b; if (is.null(dim(b_matrix))) b_matrix<-matrix(b_matrix); part2beta_4lm_help(b = b_matrix[, model_index_within_group(i)], S = SS[, model_group(i)], X = X, y = y, fl=fl)})
    heights <- sapply(rev(idx), function(i) mm[[model_group(i)]]$heights[model_index_within_group_inverted(i)])

    rownames(be) <- colnames(x.full)
    if(length(ord) > 0){
                  ind1 <- c(1:p)
                  ind1[-ord] = 1:(p - length(ord))
                  ind1[ord] = (p - length(ord) + 1):p
                  be = be[ind1,]
   }
   fit <- list(beta = be, df = length(idx):1, rss = rss[cbind(idx, ind[idx])], n = n, levels.listed = n.levels.listed, heights = heights, arguments = list(family = "gaussian", clust.method = clust.method, nlambda = nlambda, maxp = maxp), interc = TRUE)
   class(fit) = "DMR"
   return(fit)
}
