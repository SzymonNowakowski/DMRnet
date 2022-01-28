cv_compute_model<-function(model_function, Xtr, ytr, Xte, yte, real_n, lambda.full, ...) {

  #### remove from train and test columns causing data singularity
  ssd <- apply(Xtr, 2, function(x) length(unique(x)))
  singular<-which(ssd == 1)
  if (length(singular)>0) {
    Xte <- Xte[,-singular, drop=FALSE]
    Xtr <- Xtr[,-singular, drop=FALSE]
  }
  if (ncol(Xtr) == 0) {
    stop("Unable to perform cross validation. No columns in training set have any variablity in one of the folds")
  }

  model <- model_function(Xtr, ytr, ..., lambda = lambda.full)

  ###SzN remove from test the data with factors not present in training
  nn <- sapply(1:ncol(Xte), function(i) class(Xte[,i]))
  factor_columns <- which(nn == "factor")
  n.factors <- length(factor_columns)
  if (n.factors > 0)
    for (i in 1:n.factors) {

      train.levels <- model$levels.listed[[i]]

      yte<-yte[which(Xte[,factor_columns[i]] %in% train.levels), drop=FALSE]  #leaving only the test rows with levels compatible with training data
      Xte<-Xte[which(Xte[,factor_columns[i]] %in% train.levels), , drop=FALSE]

    }
  real_n <- real_n + length(yte)

  #TODO: maybe one can do better if all test data is removed in one of the folds?
  if (length(yte) == 0) {
    stop("Unable to perform cross validation. Empty test set in one of the folds")
  }

  return (list(model=model, Xtr=Xtr, ytr=ytr, Xte=Xte, yte=yte, real_n=real_n))
}
