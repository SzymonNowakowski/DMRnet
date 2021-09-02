

source("cv_helper.R")


cv_glamer <- function(X, y, family = "gaussian", clust.method = 'complete', o = 5, nlambda = 20, interc = TRUE, nfolds = 10, agressive=FALSE, maxp = ifelse(family == "gaussian", ceiling(length(y)/2), ceiling(length(y)/4))){
  X <- data.frame(X, check.names = TRUE, stringsAsFactors = TRUE)
  if (family == "gaussian"){
    n <- length(y)
    foldid <- sample(rep(1:nfolds,length.out=n))   #PP replaces nfolds by a simpler sample(rep()) function
    error <- list()

    real_n <- 0 #recount  of test instances

    glamer.full <- glamer_4lm(X, y, clust.method = clust.method, nlambda = nlambda, maxp = ceiling(maxp))
    lambdas.full<- glamer.full$lambda

    for (fold in 1:nfolds){
      Xte <- X[foldid == fold, ,drop = FALSE]
      yte <- y[foldid == fold]
      Xtr <- X[foldid != fold, ,drop = FALSE]
      ytr <- y[foldid != fold]

      helper<- cv_helper(Xtr, ytr, Xte, yte, real_n, agressive)
      Xtr<-helper$Xtr
      ytr<-helper$ytr
      Xte<-helper$Xte
      yte<-helper$yte

      glamer <- glamer_4lm(Xtr, ytr, clust.method = clust.method, lambda = lambdas.full, maxp = ceiling(maxp))
      #above: lambda - when provided - overrides nlambda

      pred <- predict(glamer, newx = as.data.frame(Xte))
      error[[fold]] <- apply(pred, 2, function(z) sum((z - yte)^2))
    }
    foldmin <- min(sapply(error, length))
    error <- sapply(1:length(error), function(i) error[[i]][(length(error[[i]]) - foldmin + 1) : length(error[[i]])])
    error <- rowSums(error)/n
    kt <- which(error == min(error))
    df.min <- glamer$df[kt[length(kt)]]
  } else{
    if (family == "binomial"){
      if (class(y) != "factor"){
        stop("Error: y should be a factor")
      }
      lev <- levels(factor(y))
      if (length(lev) != 2){
        stop("Error: factor y should have 2 levels")
      }
      n1 <- table(y)[1]
      n2 <- table(y)[2]
      foldid1 <- sample(rep(1:nfolds,length.out=n1))  #PP replaces nfolds by a simpler sample(rep()) function
      foldid2 <- sample(rep(1:nfolds,length.out=n2))  #PP replaces nfolds by a simpler sample(rep()) function
      foldid <- c()
      foldid[which(y == levels(factor(y))[1])] = foldid1
      foldid[which(y == levels(factor(y))[2])] = foldid2
      error <- list()

      real_n <- 0 #recount  of test instances

      glamer.full <- glamer_4glm(X, y, clust.method = clust.method, nlambda = nlambda, maxp = maxp)
      lambdas.full<- glamer.full$lambda

      for (fold in 1:nfolds){
        Xte <- X[foldid == fold, , drop = FALSE]
        yte <- y[foldid == fold]
        Xtr <- X[foldid != fold, , drop = FALSE]
        ytr <- y[foldid != fold]

        helper<- cv_helper(Xtr, ytr, Xte, yte, real_n, agressive)
        Xtr<-helper$Xtr
        ytr<-helper$ytr
        Xte<-helper$Xte
        yte<-helper$yte

        glamer <- glamer_4glm(Xtr, ytr, clust.method = clust.method, lambda = lambdas.full, maxp = ceiling(maxp))
        #above: lambda - when provided - overrides nlambda

        pred <- predict(glamer, newx = as.data.frame(Xte), type = "class")
        error[[fold]] <- apply(pred, 2, function(z) sum(z != yte))
      }
      foldmin <- min(sapply(error, length))
      error <- sapply(1:length(error), function(i) error[[i]][(length(error[[i]]) - foldmin + 1) : length(error[[i]])])
      error <- rowSums(error)/(n1 + n2)
      kt <- which(error == min(error))
      df.min <- glamer$df[kt[length(kt)]]
    }
    else{
      stop("Error: wrong family, should be one of gaussian, binomial")
    }
  }
  out <- list(df.min = df.min, dmr.fit = glamer.full, cvm = error, foldid = foldid)
  class(out) <- "cv.DMR"   #this is for compatibility with predict from DMR package
  return(out)
}
