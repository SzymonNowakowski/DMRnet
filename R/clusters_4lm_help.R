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
      out <- glamer_stats(c(0,betas[i1:i2]))   #appending 0 as a beta for the constrained level. Each factor has one level constrained to have beta==0
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
  ###only cut points == 0 should reference to the 1st model
  cut_point_model_reference[cut_points==0] <- 1

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
