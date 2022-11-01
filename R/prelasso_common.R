prelasso_common <- function(X, y, nn) {
#called by DMRnet4lm, DMRnet4glm,
  #        glamer_4lm, glamer_4glm
  #        DMR4lm, DMR4glm
  out <- prelasso_cont_columns(X, y)
  n <-             out$n
  nn <-            out$nn
  p.x <-           out$p.x
  X <-             out$X

  if(is.null(colnames(X))) colnames(X) <- paste("x", 1:ncol(X), sep = "")
  names(nn) <- colnames(X)

  if(sum(nn != "numeric" & nn != "factor" ) > 0){
    stop("Error: wrong data type, columns should be one of types: integer, factor, numeric")
  }

  factor_columns <- which(nn == "factor")
  n.factors <- length(factor_columns)
  n.levels <- c()
  levels.listed<-c()
  p.fac <- 0
  if (n.factors > 0){
    X[,factor_columns]<-lapply(1:n.factors, function(i) factor(X[,factor_columns[i]]))   #recalculate factors to minimal possible set
    levels.listed<-lapply(1:n.factors, function(i) levels(X[,factor_columns[i]]))
    n.levels <- sapply(1:n.factors, function(i) length(levels.listed[[i]]))
    p.fac <- sum(n.levels - 1)   #sum of factor dimensions. Again, for a factor with k levels it is taking (k-1) as a summand
  }

  cont <- which(nn == "numeric")
  n.cont <- length(cont)
  names.cont <- names(nn)[cont]

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

  X <- X[, order(nn), drop = FALSE]  # for DMR it is important that it happens early on, before factor_columns is computed
                                     # OR: compute factor_columsn once more
  nn <- sort(nn)       # I am assuming the sort IS STABLE
  factor_columns <- which(nn == "factor")       #recalculated for the sake of DMR

  x.full <- stats::model.matrix(y~., data = data.frame(y=y, X, check.names = TRUE))
  p <- ncol(x.full)   # sum over all dimensions of X, e.g. for a factor with k levels it is taking (k-1) as a summand
  fl <- c(n.levels, rep(2, n.cont))
  x.full_normalized <- apply(x.full, 2, function(x) sqrt(n/sum(x^2))*x)

  groups <- rep(1:p.x, fl-1)

  return(list(X=X, factor_columns=factor_columns, n=n, n.levels=n.levels, n.factors=n.factors, n.cont=n.cont, names.cont=names.cont, levels.listed=levels.listed, fl=fl, x.full=x.full, x.full_normalized=x.full_normalized, p=p, p.x=p.x, p.fac=p.fac, ord=ord, groups=groups))
}
