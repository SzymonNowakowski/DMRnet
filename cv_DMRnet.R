library(glmnet)

cv_helper<-function(Xtr, ytr, Xte, yte, real_n) {


  ####SzN remove from train and test columns causing data singularity
  #removing columns with only one level:
  singular_factors<-which(sapply(sapply(Xtr, levels), length)==1)  #for continous columns length is 0
  if (length(singular_factors)>0) {
    Xte <- Xte[,-singular_factors]
    Xtr <- Xtr[,-singular_factors]
    cat("removed", length(singular_factors), "columns due to singular factors in training set\n")
  }

  ###SzN remove from test the data with factors not present in training
  nn <- sapply(1:ncol(Xte), function(i) class(Xte[,i]))
  faki <- which(nn == "factor")
  n.factors <- length(faki)
  if (n.factors > 0)
    for (i in 1:n.factors) {
      Xtr[,faki[i]] <- factor(Xtr[,faki[i]])  #may be removed in next version of package because
      train.levels <- levels(Xtr[,faki[i]])   #the same info is in dmr$levels.listed[[i]]
      #but this has the advantage that is also compatible with the old package version
      yte<-yte[which(Xte[,faki[i]] %in% train.levels)]
      Xte<-Xte[which(Xte[,faki[i]] %in% train.levels),]
      Xte[,faki[i]]<-factor(Xte[,faki[i]], levels=train.levels)  #recalculate the factors in new test set, may be removed in next version of package because the same happens in predict(..)
    }
  real_n <- real_n + length(yte)
  #TODO: what to do if all test data is removed?



  #preparation to detect data dependency
  Xtr.make <- makeX(Xtr)
  prev_pos <- 0
  first_identified_yet = FALSE
  for (i in 1:ncol(Xtr))
    if (i %in% faki) {  #removing columns from the last level (but not for the first original column) they are linearly dependant
      # cat(i, prev_pos, length(levels(insurance.train.10percent.x[,i])), "\n")
      if (first_identified_yet) {
        Xtr.make <- Xtr.make[,-(prev_pos+length(levels(Xtr[,i])))]
        prev_pos <- prev_pos+length(levels(Xtr[,i])) - 1
      } else {
        prev_pos <- prev_pos+length(levels(Xtr[,i]))
        first_identified_yet = TRUE
      }
    } else prev_pos<-prev_pos+1

  QR<- qr(Xtr.make)

  if (QR$rank < ncol(Xtr.make)) {  #singular
    reverse_lookup<-rep(0, ncol(Xtr.make))
    pos<-1
    first_identified_yet = FALSE
    for (i in 1:ncol(Xtr))
      if (i %in% faki) {
        if (first_identified_yet) {
          reverse_lookup[pos:(pos+length(levels(Xtr[,i]))-2)]<-i  #there are levels-1 columns corresponding to the first original column
          pos<-pos+length(levels(Xtr[,i]))-1
        } else {
          reverse_lookup[pos:(pos+length(levels(Xtr[,i]))-1)]<-i  #there are levels-1 columns corresponding to each original column other than the first
          pos<-pos+length(levels(Xtr[,i]))
          first_identified_yet = TRUE
        }
      } else {
        reverse_lookup[pos]<-i
        pos<-pos+1
      }

    #removal of columns for pivot positions larger than rank
    remove_us<-reverse_lookup[QR$pivot[(QR$rank+1):length(QR$pivot)]]
    Xtr <- Xtr[,-remove_us]
    Xte <- Xte[,-remove_us]
    cat("removed", length(unique(remove_us)), "columns due to data linear dependency\n")
  }
  return (list(Xtr=Xtr, ytr=ytr, Xte=Xte, yte=yte, real_n=real_n))
}

cv_DMRnet <- function(X, y, family = "gaussian", clust.method = 'complete', o = 5, nlambda = 20, lam = 10^(-7), interc = TRUE, nfolds = 10, method = "DMRnet", maxp = ifelse(family == "gaussian", ceiling(length(y)/2), ceiling(length(y)/4))){


  X <- data.frame(X, check.names = TRUE, stringsAsFactors = TRUE)
  if (family == "gaussian"){
    n <- length(y)
    real_n <- 0 #recount  of test instances
    #PP new code foldid <- cvfolds(n, nfolds)
    foldid <- sample(rep(1:nfolds,length.out=n))   #PP replaces nfolds by a simpler sample(rep()) function
    #PP new code error <- list()
    err <- list(); rss <- list(); #md <- list()

    for (fold in 1:nfolds){
      cat("gaussian fold:", fold, "\n")
      if (fold==2) {
        xx<-4
      }
      Xte <- X[foldid == fold, ,drop = FALSE]
      yte <- y[foldid == fold]
      Xtr <- X[foldid != fold, ,drop = FALSE]
      ytr <- y[foldid != fold]

      helper<- cv_helper(Xtr, ytr, Xte, yte, real_n)
      Xtr<-helper$Xtr
      ytr<-helper$ytr
      Xte<-helper$Xte
      yte<-helper$yte
      real_n<-helper$real_n

      if (method == "GLAF") {
        dmr <- glaf_4lm(Xtr, ytr, clust.method = clust.method, nlambda = nlambda, maxp = ceiling(maxp))
        cat("GLAF calculated\n")

      } else dmr <- DMRnet(Xtr, ytr, family = "gaussian", clust.method = clust.method, o = o, nlambda = nlambda, interc = interc, maxp = ceiling(maxp))
      #PP new code
      rss[[fold]] <- dmr$rss

      pred <- predict(dmr, newx = as.data.frame(Xte))
      #PP new code error[[fold]] <- apply(pred, 2, function(z) sum((z - yte)^2))
      err[[fold]] <- apply(pred, 2, function(z) mean((z - yte)^2))

    }

    #PP new code foldmin <- min(sapply(error, length))
    #            error <- sapply(1:length(error), function(i) error[[i]][(length(error[[i]]) - foldmin + 1) : length(error[[i]])])
    #            error <- rowSums(error)/n
    len_err <- sapply(err, length)
    foldmin <- min(len_err)
    ERR <- sapply(1:nfolds, function(i) err[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
    #err <- rowMeans(ERR); kt <- which(err == min(err)); df.min <- dmr$df[kt[length(kt)]]; plot(err, type="o")

    #PP rename dmr.fit
    if (method == "GLAF") {
      dmr.full <- glaf_4lm(X, y, clust.method = clust.method, nlambda = nlambda, maxp = ceiling(maxp))
      cat("full GLAF calculated\n")
    } else dmr.full <- DMRnet(X, y, family = "gaussian", clust.method = clust.method, o = o, nlambda = nlambda, interc = interc, maxp = ceiling(maxp))

    #PP new code kt <- which(error == min(error))
    #            df.min <- dmr$df[kt[length(kt)]]
    p1 <- dmr.full$df[1]
    s2 <- dmr.full$rss[1]/(n-p1)
    Const <- exp(seq(log(2/50),log(2*50), length=80))
    laGIC <- Const*log(p1)*s2
    RSS <- sapply(1:nfolds, function(i) rss[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
    #MD <- sapply(1:nfolds, function(i)  md[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
    IND <- apply( RSS, 2, function(r) sapply( laGIC, function(la) which.min(r+la*length(r):1) ) )
    errGIC <- apply( IND, 1, function(ind) mean(ERR[cbind(ind,1:nfolds)]) )
    #mdGIC  <- apply( IND, 1, function(ind) mean(MD[cbind(ind,1:10)]) )
    #plot(mdGIC[length(laGIC):1],errGIC[length(laGIC):1]/s2, xlab="MD", ylab="PE", type="o")

    r <- dmr.full$rss
    kt <- which(errGIC == min(errGIC))
    indGIC <- kt[length(kt)]
    gic.full <- (r+laGIC[indGIC]*length(r):1)/(real_n*s2)
    #plot(gic.full[length(gic.full):1])
    indMod <- which.min(gic.full)
    df.min <- dmr.full$df[indMod]



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
      real_n <- 0 #recount  of test instances

      foldid1 <- sample(rep(1:nfolds,length.out=n1))  #PP replaces nfolds by a simpler sample(rep()) function
      foldid2 <- sample(rep(1:nfolds,length.out=n2))  #PP replaces nfolds by a simpler sample(rep()) function
      foldid <- c()
      foldid[which(y == levels(factor(y))[1])] = foldid1
      foldid[which(y == levels(factor(y))[2])] = foldid2
      #PP new code error <- list()
      err <- list(); loglik <- list(); #md <- list()
      for (fold in 1:nfolds) {
        cat("binomial fold:", fold, "\n")
        Xte <- X[foldid == fold, ,drop = FALSE]
        yte <- y[foldid == fold]
        Xtr <- X[foldid != fold, ,drop = FALSE]
        ytr <- y[foldid != fold]

        helper<- cv_helper(Xtr, ytr, Xte, yte, real_n)
        Xtr<-helper$Xtr
        ytr<-helper$ytr
        Xte<-helper$Xte
        yte<-helper$yte
        real_n<-helper$real_n

        if (method == "GLAF") {
          dmr <- glaf_4glm(Xtr, ytr, clust.method = clust.method, nlambda = nlambda, maxp = ceiling(maxp))
          cat("GLAF calculated\n")

        } else dmr <- DMRnet(Xtr, ytr, family = "binomial", clust.method = clust.method, o = o, nlambda = nlambda, lam = lam, interc = interc, maxp = maxp)
        #SzN new code based on PP new code
        loglik[[fold]] <- -2*dmr$loglik

        pred <- predict(dmr, newx = as.data.frame(Xte), type = "class")
        #SzN new code based on PP new code error[[fold]] <- apply(pred, 2, function(z) sum(z != yte))
        err[[fold]] <- apply(pred, 2, function(z) mean(z != yte))

      }

      #PP new code foldmin <- min(sapply(error, length))
      #            error <- sapply(1:length(error), function(i) error[[i]][(length(error[[i]]) - foldmin + 1) : length(error[[i]])])
      #            error <- rowSums(error)/(n1+n2)
      len_err <- sapply(err, length)
      foldmin <- min(len_err)
      ERR <- sapply(1:nfolds, function(i) err[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
      #err <- rowMeans(ERR); kt <- which(err == min(err)); df.min <- dmr$df[kt[length(kt)]]; plot(err, type="o")


      #PP rename dmr.fit
      if (method == "GLAF") {
        dmr.full <- glaf_4glm(X, y, clust.method = clust.method, nlambda = nlambda, lam = lam, maxp = maxp)
        cat("full GLAF calculated\n")

      } else dmr.full <- DMRnet(X, y, family = "binomial", clust.method = clust.method, o = o, nlambda = nlambda, lam = lam, interc = interc, maxp = maxp)

      #SzN new code based on PP's new code kt <- which(error == min(error))
      #            df.min <- dmr$df[kt[length(kt)]]
      p1 <- dmr.full$df[1]
      Const <- exp(seq(log(2/50),log(2*50), length=80))
      laGIC <- Const*log(p1)
      LOGLIK <- sapply(1:nfolds, function(i) loglik[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
      #MD <- sapply(1:nfolds, function(i)  md[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
      IND <- apply( LOGLIK, 2, function(ll) sapply( laGIC, function(la) which.min(ll+la*length(ll):1) ) )
      errGIC <- apply( IND, 1, function(ind) mean(ERR[cbind(ind,1:nfolds)]) )
      #mdGIC  <- apply( IND, 1, function(ind) mean(MD[cbind(ind,1:10)]) )
      #plot(mdGIC[length(laGIC):1],errGIC[length(laGIC):1]/s2, xlab="MD", ylab="PE", type="o")

      ll <- -2*dmr.full$loglik
      kt <- which(errGIC == min(errGIC))
      indGIC <- kt[length(kt)]
      gic.full <- (ll+laGIC[indGIC]*length(ll):1)/real_n
      #plot(gic.full[length(gic.full):1])
      indMod <- which.min(gic.full)
      df.min <- dmr.full$df[indMod]

    }
    else{
      stop("Error: wrong family, should be one of gaussian, binomial")
    }
  }
  #PP: out <- list(df.min = df.min, dmr.fit = dmr.fit, cvm = error, foldid = foldid)
  out <- list(df.min = df.min, dmr.fit = dmr.full, cvm = errGIC, foldid = foldid)
  class(out) <- "cv.DMR"
  cat("cv.out\n")
  return(out)
}
