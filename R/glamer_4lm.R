glamer_4lm <- function(X, y, cut_points = NULL, clust.method = "complete", lambda = NULL, nlambda = 100, maxp = ceiling(length(y)/2)){

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

  if (is.null(lambda)) {
    user.lambdas<-substitute()    #make user.lambdas - paradoxically - not present in a call to grpreg
  } else {
    nlambda <- length(lambda)   #override this parameter
    user.lambdas <- lambda
  }

  mL <- grpreg::grpreg(x.full[,-1], y, group=rep(1:p.x, fl-1) , penalty = "grLasso", family ="gaussian", nlambda = nlambda, lambda = user.lambdas)
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

  shift <- function(i) {sum(rss[, i] == Inf)}

  if (length(cut_points)>0) {
    cut_point_model_reference <- sapply(1:length(mm), function(i) mm[[i]]$cut_point_model_reference)
    mask <- matrix(Inf, nrow=nrow(rss), ncol=ncol(rss))   #it will mask by Inf all models not back-referenced in cut_points
    model_cut_point_reference_max <- model_cut_point_reference_min <- matrix(0, nrow=nrow(rss), ncol=ncol(rss))
    for (i in 1:length(mm)) {
      mask[shift(i) + na.omit(cut_point_model_reference[,i]), i] <- 0

      for (model in 1:length(mm[[i]]$rss)) {
        model_cut_point_reference_min[shift(i) + model, i] <- which(cut_point_model_reference[,i]==model)[1]   #first occurence of TRUE
        model_cut_point_reference_max[shift(i) + model, i] <- length(cut_point_model_reference[,i]) - which(rev(cut_point_model_reference[,i])==model)[1] +1   #last occurence of TRUE
      }
      # cut_point_model_reference, for each cutpoint, it lists models the cutpoint references
      # on the other hand model_cut_point_reference is a reverse relation:
      # model_cut_point_reference, for each model, it lists the cutpoints that originate the models (min and max such cutpoint)
      # both 0 and NA mean "no related cutpoint"
    }
    rss_with_mask <- rss + mask
  } else
    rss_with_mask <- rss
  ind <- apply(rss_with_mask, 1, which.min)  #in each row, which is the index of a model minimizing rss
  indInf <- apply(rss_with_mask,1,min) == Inf



  maxi <- min(p, maxp)
  if (length(ind) > maxi){
    idx <- (length(ind) - maxi):length(ind)   #but real model sizes are still length(idx):1
  } else {
    idx <- 1:length(ind)
  }

  #smallest models are last
  model_group <- function(i) {ind[i]}
  model_index_within_group <- function(i) {i - shift(model_group(i))}
  #model_index_within_group_inverted <- function(i) {length(mm[[model_group(i)]]$heights)-model_index_within_group(i)+1}


  be <- sapply(idx, function(i) {
    b_matrix<-mm[[model_group(i)]]$b;
    if (is.null(dim(b_matrix))) {
      b_matrix<-matrix(b_matrix);    #TODO: verify if this shouldn't be b_matrix<-t(as.matrix(b_matrix))   as with other VERTICAL matrices that degenerated to HORIZONTAL vectors
    }
    part2beta_help(b = b_matrix[, model_index_within_group(i)], S = SS[, model_group(i)], X = X, y = y, fl=fl, valid=!indInf[i])
  })

  heights <- sapply(idx, function(i) mm[[model_group(i)]]$heights[model_index_within_group(i)])
  #heights from full model #1 with height == 0 to the last 1-element model with height > 0

  if (length(cut_points>0)) {
    cut_point_generated_models_min <- sapply(idx, function(i) {ifelse(indInf[i],0,model_cut_point_reference_min[i, model_group(i)])})
    cut_point_generated_models_max <- sapply(idx, function(i) {ifelse(indInf[i],0,model_cut_point_reference_max[i, model_group(i)])})
    rss_vec<-sapply(idx, function(i) {ifelse(indInf[i],Inf,rss[i, model_group(i)])})
  }
  else {
    cut_point_generated_models_min <- cut_point_generated_models_max <- rep(0, length(idx))
    rss_vec<-rss[cbind(idx, ind[idx])]
  }

  rownames(be) <- colnames(x.full)
  if(length(ord) > 0){
    ind1 <- c(1:p)
    ind1[-ord] = 1:(p - length(ord))
    ind1[ord] = (p - length(ord) + 1):p
    be = be[ind1,]
  }

  fit <- list(beta = be, df = length(idx):1, rss = rss_vec, n = n, levels.listed = n.levels.listed, heights = heights, lambda = mL$lambda, cut_point_generated_models = list(min = cut_point_generated_models_min, max = cut_point_generated_models_max), arguments = list(family = "gaussian", clust.method = clust.method, nlambda = nlambda, maxp = maxp), interc = TRUE)

  class(fit) = "DMR"
  return(fit)
}
