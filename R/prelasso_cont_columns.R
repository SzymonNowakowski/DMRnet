prelasso_cont_columns <- function(X, y) {
# called by SOSnet4lm, SOSnet4glm and prelasso_common
  n <- nrow(X)
  if(n != length(y)){
    stop("Error: non-conforming data: nrow(X) not equal to length(y)")
  }
  ssd <- apply(X, 2, function(x) length(unique(x)))
  if (ssd[1] == 1 & (inherits(X[,1], "numeric") | inherits(X[,1], "integer"))){ # removing the first column
    # in case it is a numeric constant in X
    # i.e. in case it is an Intercept. Other than that, constant columns are NOT allowed
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

  nn[nn == "integer"] <- "numeric"
  p.x <- ncol(X)
  return(list(X=X, n=n, nn=nn, p.x=p.x))
}
