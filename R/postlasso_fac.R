postlasso_fac <- function(bb, groups) {

  bb <- bb[-1, ,drop=FALSE]   #fac must be computed on bb without intercept

  fac <- apply(bb, 2, function(x) tapply(x, factor(groups), function(z) sum(z^2)*sqrt(length(z))))

#fac is a matrix with normalized beta statistics relating to GROUPS of variables (rows) respective to non-duplicated lambdas (colums)
#nrow = #GROUPS of variables
#ncol = #active lambdas
  if(is.null(dim(fac))){  #in case of a single k-level factor matrix in X, there is only one group and fac would be reduced to a vector. This line here helps to convert it back to a matrix
  #unfortunately a symmetric situation is still possible although grpreg does NOT accept a single lambda value nor nlambda=1, nlambda must be at least two
  #but lamdas get removed in case of having betas number too large, too little (==0) or duplicated
  #still, even if there is only one column left (only-intercept-lambda related) fac is a one column matrix, which is OK

    fac <- t(as.matrix(fac))
  }
  return(fac)
}
