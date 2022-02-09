glamer_4lm <- function(X, y, clust.method, nlambda, lam, maxp, lambda){

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
  factor_columns <- which(nn == "factor")
  n.factors <- length(factor_columns)
  n.levels <- c()
  if (n.factors > 0){
    X[,factor_columns]<-lapply(1:n.factors, function(i) factor(X[,factor_columns[i]]))   #recalculate factors
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

  groups <- rep(1:p.x, fl-1)
  #mL <- grpreg::grpreg(x.full[,-1], y, group=rep(1:p.x, fl-1) , penalty = "grLasso", family ="gaussian", nlambda = nlambda, lambda = user.lambda)
  #to keep the x.full as a matrix after the intercept column is removed in case of a single 2-level factor column in X, drop=FALSE must be provided in grpreg call below in x.full[,-1, drop=FALSE]
  #otherwise x.full[,-1] would be reduced to a vector and grpreg would crash
  mL <- grpreg::grpreg(x.full[,-1, drop=FALSE], y, group=groups , penalty = "grLasso", family ="gaussian", nlambda = nlambda, lambda = user.lambda)
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

  fac <- apply(bb[-1,ii == FALSE, drop = FALSE], 2, function(x) tapply(x, factor(groups), function(z) sum(z^2)*sqrt(length(z))))

  if(is.null(dim(fac))){#in case of a single k-level factor matrix in X, there is only one group and fac would be reduced to a vector. This line here helps to convert it back to a matrix
    #by the way, a symmetric situation is not possible as grpreg does NOT accept a single lambda value nor nlambda=1, nlambda must be at least two
    fac <- t(as.matrix(fac))
  }

  bb[1+ which(groups %in% which(fac>0)),]<-bb[1+ which(groups %in% which(fac>0)),]+lam
            #regularizing those betas that belong to groups > 0 to make DOUBLE SURE THEY ARE STRICTLY >0
            #(there have been cases of grpreg not observing the group constraint - vide hard_case_DMRnet_promoter test file in testing_branch)


  SS <- sapply(1:sum(ii == FALSE), function(i) ifelse(fac[, i] > 0, 1, 0))
  if(is.null(dim(SS))){   #for a single variable (a single k-level factor) SS is a vector
    SS <- t(as.matrix(SS))  #change it into a HORIZONTAL matrix
  }

  #if (p >= n) SS = SS[,-1]  #in original DMRnet code, this line serves to eliminate the first FULL model if p >=n, because it would cause problems in DMR later. However, this is not the case here in GLAMER
                            #because the SS matrix is built in a different way and doesn't include the full model as the first column.



  bb<-bb[,ii==FALSE, drop=FALSE]   #betas but for active lambdas only

  mm <- lapply(1:ncol(SS), function(i) glamer_4lm_help(SS[,i], bb[,i], mL, X, y, fl, clust.method, lam = lam))
  maxl <- max(sapply(1:length(mm), function(i) length(mm[[i]]$rss)))
  rss <- sapply(1:length(mm), function(i) c(rep(Inf, maxl - length(mm[[i]]$rss)), mm[[i]]$rss))

  ind <- apply(rss, 1, which.min)  #in each row, which is the index of a model minimizing rss

  maxi <- min(p, maxp)
  if (length(ind) > maxi){
    idx <- (length(ind) - maxi):length(ind)   #but real model sizes are still length(idx):1
  } else {
    idx <- 1:length(ind)
  }

  #smallest models are last
  shift<- function(i) {sum(rss[, i]==Inf)}
  model_group <- function(i) {ind[i]}
  model_index_within_group<- function(i) {i-shift(model_group(i))}

  be <- sapply(idx, function(i) {
    b_matrix<-mm[[model_group(i)]]$b;
    if (is.null(dim(b_matrix))) {
      b_matrix<-matrix(b_matrix);    #note this shouldn't be b_matrix<-t(as.matrix(b_matrix)) :   with other matrices that degenerated to HORIZONTAL vectors we want to have HORIZONTAL matrices
      #in this case, however, we want a VERTICAL matrix
      #in b_matrix[,model_index_within_group(i)] we take columns of b_matrix if b_matrix is a legitimate matrix
      #when it is degenerate (a vector), we want that taking the full first column (as model_index_within_group(i) == 1 in those cases) takes this whole vector
    }
    part2beta_help(b = b_matrix[, model_index_within_group(i)], S = SS[, model_group(i)], X = X, y = y, fl=fl)
  })

  rownames(be) <- colnames(x.full)
  if(length(ord) > 0){
    ind1 <- c(1:p)
    ind1[-ord] = 1:(p - length(ord))
    ind1[ord] = (p - length(ord) + 1):p
    be = be[ind1,]
  }

  fit <- list(beta = be, df = length(idx):1, rss = rss[cbind(idx, ind[idx])], n = n, levels.listed = levels.listed, lambda = mL$lambda, arguments = list(family = "gaussian", clust.method = clust.method, nlambda = nlambda, maxp = maxp), interc = TRUE)

  class(fit) = "DMR"
  return(fit)
}
