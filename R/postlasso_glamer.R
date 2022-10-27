postlasso_glamer <- function(bb, lam, fac, groups) {


  #first, note that sum(ii==FALSE) is a number of predictor sets and it may be smaller than nlambda because of (1), (2), (3)
  SS <- sapply(1:ncol(bb), function(i) ifelse(fac[, i] > 0, 1, 0))
  #nrow = #predictors
  #ncol = #active lambdas
  if(is.null(dim(SS))){   #for a single variable (a single k-level factor) SS is a vector
    SS <- t(as.matrix(SS))  #change it into a matrix with one row - HORIZONTAL matrix is the effect
  }

  #if (p >= n) SS <- SS[,-1]   #in original DMRnet code, this line serves to eliminate the first FULL model if p >=n, because it would cause problems in DMR later. However, this is not the case here in GLAMER
  #because the SS matrix is built in a different way and doesn't include the full model as the first column.

  bb<-rbind(bb[1, , drop=FALSE], sapply(1:ncol(bb), function(c) sapply(1:(nrow(bb)-1), function(r) bb[1+r,c]+lam*(fac[groups[r], c]>0))))
  #betas but for active lambdas only
  #regularizing those betas that belong to groups > 0 to make DOUBLE SURE THEY ARE STRICTLY >0
  #(there have been cases of grpreg not observing the group constraint - vide hard_case_DMRnet_promoter test file in testing_branch)

  return(list(bb=bb, SS=SS))
}
