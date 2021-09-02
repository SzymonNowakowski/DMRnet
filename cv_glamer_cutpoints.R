library(glmnet) #for makeX function

source("cv_helper.R")

cv_glamer_cutpoints <- function(X, y, family = "gaussian", clust.method = 'complete', o = 5, nlambda = 20, lam = 10^(-7), interc = TRUE, nfolds = 10, agressive = FALSE, maxp = ifelse(family == "gaussian", ceiling(length(y)/2), ceiling(length(y)/4))){


  X <- data.frame(X, check.names = TRUE, stringsAsFactors = TRUE)
  if (family == "gaussian"){
    n <- length(y)
    real_n <- 0 #recount  of test instances
    #PP new code foldid <- cvfolds(n, nfolds)
    foldid <- sample(rep(1:nfolds,length.out=n))   #PP replaces nfolds by a simpler sample(rep()) function
    #PP new code error <- list()
    #err <- list(); #rss <- list(); #md <- list()

    glamer.full <- glamer_4lm(X, y, clust.method = clust.method, nlambda = nlambda, maxp = ceiling(maxp))
    lambdas.full<- glamer.full$lambda
    heights<- sort(glamer.full$heights)

    first_height <- heights[2] / 10
    heights_shifted <- c(first_height, heights[2:(length(heights)-1)])
    cut_points <- c(0, sqrt(heights[2:length(heights)] * heights_shifted), heights[length(heights)])   #geometric mean, equivalent to arithmetic mean on a logarithmic scale
    error_per_cutpoint <- lapply(1:length(cut_points), function(i) {numeric(0)})   #initiating it with empty vectors

    for (fold in 1:nfolds){
      cat("gaussian fold:", fold, "\n")

      Xte <- X[foldid == fold, ,drop = FALSE]
      yte <- y[foldid == fold]
      Xtr <- X[foldid != fold, ,drop = FALSE]
      ytr <- y[foldid != fold]

      helper<- cv_helper(Xtr, ytr, Xte, yte, real_n, agressive)
      Xtr<-helper$Xtr
      ytr<-helper$ytr
      Xte<-helper$Xte
      yte<-helper$yte
      real_n<-helper$real_n


      glamer <- glamer_4lm(Xtr, ytr, cut_points, clust.method = clust.method, lambda = lambdas.full, maxp = ceiling(maxp))

      pred <- predict(glamer, newx = as.data.frame(Xte))
      #PP new code error[[fold]] <- apply(pred, 2, function(z) sum((z - yte)^2))
      err<- apply(pred, 2, function(z) mean((z - yte)^2))

      for (m_size in 1:length(glamer$rss))
        if (glamer$rss[m_size]<Inf)   #there was a cutpoint associated with that model size
          for (cutpoint in glamer$cut_point_generated_models$min[m_size] : glamer$cut_point_generated_models$max[m_size])
            error_per_cutpoint[[cutpoint]] <- c(error_per_cutpoint[[cutpoint]], err[m_size]/length(yte))

    }


    mean_error_per_cutpoint <- sapply(error_per_cutpoint, mean)
    kt <- which(mean_error_per_cutpoint <= min(na.omit(mean_error_per_cutpoint)) + sd(na.omit(mean_error_per_cutpoint)))
    minimal_higher_height = min(glamer.full$heights[glamer.full$heights >= cut_points[kt[length(kt)]]])
    kt <- which(glamer.full$heights== minimal_higher_height)
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

        helper<- cv_helper(Xtr, ytr, Xte, yte, real_n, agressive)
        Xtr<-helper$Xtr
        ytr<-helper$ytr
        Xte<-helper$Xte
        yte<-helper$yte
        real_n<-helper$real_n

        if (method == "GLAMER") {
          dmr <- glamer_4glm(Xtr, ytr, clust.method = clust.method, nlambda = nlambda, maxp = ceiling(maxp))
          cat("GLAMER calculated\n")

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
      if (method == "GLAMER") {
        glamer.full <- glamer_4glm(X, y, clust.method = clust.method, nlambda = nlambda, lam = lam, maxp = maxp)
        cat("full GLAMER calculated\n")

      } else glamer.full <- DMRnet(X, y, family = "binomial", clust.method = clust.method, o = o, nlambda = nlambda, lam = lam, interc = interc, maxp = maxp)

      #SzN new code based on PP's new code kt <- which(error == min(error))
      #            df.min <- dmr$df[kt[length(kt)]]
      p1 <- glamer.full$df[1]
      Const <- exp(seq(log(2/50),log(2*50), length=80))
      laGIC <- Const*log(p1)
      LOGLIK <- sapply(1:nfolds, function(i) loglik[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
      #MD <- sapply(1:nfolds, function(i)  md[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
      IND <- apply( LOGLIK, 2, function(ll) sapply( laGIC, function(la) which.min(ll+la*length(ll):1) ) )
      errGIC <- apply( IND, 1, function(ind) mean(ERR[cbind(ind,1:nfolds)]) )
      #mdGIC  <- apply( IND, 1, function(ind) mean(MD[cbind(ind,1:10)]) )
      #plot(mdGIC[length(laGIC):1],errGIC[length(laGIC):1]/s2, xlab="MD", ylab="PE", type="o")

      ll <- -2*glamer.full$loglik
      kt <- which(errGIC == min(errGIC))
      indGIC <- kt[length(kt)]

      gic.full <- (ll+laGIC[indGIC]*length(ll):1)/real_n
      #plot(gic.full[length(gic.full):1])
      # Plateau-resistant CV:
      if (plateau_resistant_CV) {
        indMods <- which(gic.full <= min(gic.full) + sd(gic.full))
        indMod <- indMods[length(indMods)]
      } else
        indMod <- which.min(gic.full)

      df.min <- glamer.full$df[indMod]

    }
    else{
      stop("Error: wrong family, should be one of gaussian, binomial")
    }
  }
  #PP: out <- list(df.min = df.min, dmr.fit = dmr.fit, cvm = error, foldid = foldid)
  out <- list(df.min = df.min, dmr.fit = glamer.full, cvm = error_per_cutpoint, foldid = foldid)
  class(out) <- "cv.DMR"
  cat("cv.out\n")
  return(out)
}
