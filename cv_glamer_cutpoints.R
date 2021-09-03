library(glmnet) #for makeX function

source("cv_helper.R")


cut_points_from_heights <- function(heights) {
  #heights should be >0 and may have 0.0 as one of the values
  #resulting cut_points,, they always have 0 (full model) as one of the cut-points.
  heights <- sort(unique(heights))
  if (heights[1] == 0 & length(heights) >= 3) {
    #geometric mean is equivalent to arithmetic mean on a logarithmic scale
    heights_shifted <- heights[2:(length(heights)-1)]
    geometric_means <- sqrt(heights[3:length(heights)] * heights_shifted)
  } else if (heights[1] > 0 & length(heights) >= 2) {
    heights_shifted <- heights[1:(length(heights)-1)]
    geometric_means <- sqrt(heights[2:length(heights)] * heights_shifted)
  } else
    geometric_means <- numeric(0)

  return(sort(unique(c(0, heights, geometric_means))))

}

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
    cut_points<-cut_points_from_heights(glamer.full$heights)

    error_per_cutpoint <- lapply(1:length(cut_points), function(i) {numeric(0)})   #initiating it with empty vectors

    for (fold in 1:nfolds){
      cat("gaussian fold:", fold, "\n")

      Xte <- X[foldid == fold, ,drop = FALSE]
      yte <- y[foldid == fold]
      Xtr <- X[foldid != fold, ,drop = FALSE]
      ytr <- y[foldid != fold]

      helper <- cv_helper(Xtr, ytr, Xte, yte, real_n, agressive)
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
            error_per_cutpoint[[cutpoint]] <- c(error_per_cutpoint[[cutpoint]], err[m_size])

    }




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
      #err <- list(); loglik <- list(); #md <- list()

      glamer.full <- glamer_4glm(X, y, clust.method = clust.method, nlambda = nlambda, maxp = ceiling(maxp))

      lambdas.full<- glamer.full$lambda
      cut_points<-cut_points_from_heights(glamer.full$heights)

      error_per_cutpoint <- lapply(1:length(cut_points), function(i) {numeric(0)})   #initiating it with empty vectors

      for (fold in 1:nfolds) {
        cat("binomial fold:", fold, "\n")
        Xte <- X[foldid == fold, ,drop = FALSE]
        yte <- y[foldid == fold]
        Xtr <- X[foldid != fold, ,drop = FALSE]
        ytr <- y[foldid != fold]

        helper <- cv_helper(Xtr, ytr, Xte, yte, real_n, agressive)
        Xtr<-helper$Xtr
        ytr<-helper$ytr
        Xte<-helper$Xte
        yte<-helper$yte
        real_n<-helper$real_n

        glamer <- glamer_4glm(Xtr, ytr, cut_points, clust.method = clust.method, lambda = lambdas.full, maxp = ceiling(maxp))

        pred <- predict(glamer, newx = as.data.frame(Xte), type = "class")

        err <- apply(pred, 2, function(z) mean(z != yte))

        for (m_size in 1:length(glamer$loglik))
          if (glamer$loglik[m_size]<Inf)   #there was a cutpoint associated with that model size
            for (cutpoint in glamer$cut_point_generated_models$min[m_size] : glamer$cut_point_generated_models$max[m_size])
              error_per_cutpoint[[cutpoint]] <- c(error_per_cutpoint[[cutpoint]], err[m_size])

      }



    }
    else{
      stop("Error: wrong family, should be one of gaussian, binomial")
    }
  }

  mean_error_per_cutpoint <- sapply(error_per_cutpoint, mean)
  kt <- which(mean_error_per_cutpoint <= min(na.omit(mean_error_per_cutpoint)) + sd(na.omit(mean_error_per_cutpoint)))
  minimal_higher_height = min(glamer.full$heights[glamer.full$heights >= cut_points[kt[length(kt)]]])
  kt <- which(glamer.full$heights == minimal_higher_height)
  df.min <- glamer$df[kt[length(kt)]]

  out <- list(df.min = df.min, dmr.fit = glamer.full, cvm = error_per_cutpoint, foldid = foldid)
  class(out) <- "cv.DMR"
  cat("cv.out\n")
  return(out)
}
