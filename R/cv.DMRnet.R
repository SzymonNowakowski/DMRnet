#' @title cross-validation for DMRnet
#'
#' @description Does k-fold cross-validation for DMR and returns a value for df.
#'
#' @param X Input data frame, of dimension n x p; each row is an observation vector. Columns can be numerical or integer for continuous predictors or factors for categorical predictors.
#'
#' @param y Response variable. Numerical for family="gaussian" or a factor with two levels for family="binomial". For family="binomial" the last level in alphabetical order is the target class.
#'
#' @param family Response type; one of: "gaussian", "binomial".
#'
#' @param clust.method Clustering method used for partitioning levels of factors; see function \href{https://stat.ethz.ch/R-manual/R-devel/library/stats/html/hclust.html}{hclust} in package \pkg{stats} for details.
#'
#' @param o Parameter of the group lasso screening step, described in \code{\link{DMRnet}}.
#'
#' @param nlambda Parameter of the group lasso screening step, described in \code{\link{DMRnet}}.
#'
#' @param lam Value of parameter lambda controling the amount of penalization in rigde regression. Used only for logistic regression in order to allow for parameter estimation in linearly separable setups. Used only for numerical reasons.
#'
#' @param interc Should intercept(s) be fitted (default=TRUE) or set to zero (FALSE). If in X there are any categorical variables, interc=TRUE.
#'
#' @param nfolds Number of folds in cross-validation.
#'
#' @param maxp Maximal number of parameters of the model, smaller values result in quicker computation.
#'
#' @details cv.DMRnet algorithm does k-fold cross-validation for DMRnet. The df for the minimal estimated prediction error is returned.
#'
#' @return An object with S3 class "cv.DMR" is  returned,  which  is  a  list  with  the  ingredients  of  the  cross-validation fit.
#' \describe{
#'   \item{df.min}{df (number of parameters) for the model with minimal cross-validated error.}
#'   \item{dmr.fit}{Fitted DMR object for the full data.}
#'   \item{cvm}{The mean cross-validated error for the entire sequence of models.}
#'   \item{foldid}{The fold assignments used.}
#' }
#'
#' @seealso  \code{\link{plot.cv.DMR}} for plotting, \code{\link{coef.cv.DMR}} for extracting coefficients and \code{\link{predict.cv.DMR}} for prediction.
#'
#' @examples
#' ## cv.DMRnet for linear regression
#' set.seed(13)
#' data(miete)
#' ytr <- miete$rent[1:1500]
#' Xtr <- miete$area[1:1500]
#' Xte <- miete$area[1501:2053]
#' cv <- cv.DMRnet(Xtr, ytr)
#' print(cv)
#' plot(cv)
#' coef(cv)
#' ypr <- predict(cv, newx = Xte)
#'
#' @export cv.DMRnet

cv.DMRnet <- function(X, y, family = "gaussian", clust.method = 'complete', o = 5, nlambda = 20, lam = 10^(-7), interc = TRUE, nfolds = 10, maxp = ifelse(family == "gaussian", ceiling(length(y)/2), ceiling(length(y)/4))){


       X <- data.frame(X, check.names = TRUE, stringsAsFactors = TRUE)
       if (family == "gaussian"){
          n <- length(y)
          real_n <- 0 #recount  of test instances
          #PP new code foldid <- cvfolds(n, nfolds)
          foldid <- sample(rep(1:nfolds,length.out=n))   #PP replaces nfolds by a simpler sample(rep()) function
          #PP new code error <- list()
          err <- list(); rss <- list(); #md <- list()

          for (fold in 1:nfolds){

              Xte <- X[foldid == fold, ,drop = FALSE]
              yte <- y[foldid == fold]
              Xtr <- X[foldid != fold, ,drop = FALSE]
              ytr <- y[foldid != fold]
              dmr <- DMRnet(Xtr, ytr, family = "gaussian", clust.method = clust.method, o = o, nlambda = nlambda, interc = interc, maxp = ceiling(maxp))
              #PP new code
              rss[[fold]] <- dmr$rss

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



              pred <- predict.DMR(dmr, newx = as.data.frame(Xte))
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
          dmr.full <- DMRnet(X, y, family = "gaussian", clust.method = clust.method, o = o, nlambda = nlambda, interc = interc, maxp = ceiling(maxp))

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
              cat("fold:", fold, "\n")
              Xte <- X[foldid == fold, , drop = FALSE]
              yte <- y[foldid == fold]
              Xtr <- X[foldid != fold, , drop = FALSE]
              ytr <- y[foldid != fold]
              dmr <- DMRnet(Xtr, ytr, family = "binomial", clust.method = clust.method, o = o, nlambda = nlambda, lam = lam, interc = interc, maxp = maxp)
              #SzN new code based on PP new code
              loglik[[fold]] <- -2*dmr$loglik

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

              pred <- predict.DMR(dmr, newx = as.data.frame(Xte), type = "class")
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
          dmr.full <- DMRnet(X, y, family = "binomial", clust.method = clust.method, o = o, nlambda = nlambda, lam = lam, interc = interc, maxp = maxp)

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
