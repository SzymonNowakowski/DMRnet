SOSnet4glm <- function(X, y, o, nlambda, lam, interc, maxp, lambda){
          if (!inherits(y, "factor")){
             stop("Error: y should be a factor")
          }
          lev <- levels(y)
          if (length(lev) != 2){
             stop("Error: factor y should have 2 levels")
          }
          y <- ifelse(y == lev[2], 1, 0)
          n <- nrow(X)
          if(n != length(y)){
              stop("Error: non-conforming data: nrow(X) not equal to length(y)")
          }
          ssd <- apply(X, 2, stats::sd)
          if (ssd[1] == 0){
             X <- X[,-1, drop = FALSE]
             ssd <- ssd[-1]
          }
          if(ncol(X) == 0){
              stop("Error: X has zero columns")
          }
          if(sum(ssd == 0) > 0){
               stop("Error: X has columns with sd = 0 apart from the intercept")
          }
          nn <- sapply(1:ncol(X), function(i) class(X[,i]))
          nn[nn == "integer"] <- "numeric"
          if(sum(nn != "numeric") > 0){
              stop("Error: wrong data type, columns should be one of types: integer, numeric")
          }
          p <- ncol(X)
          Xg <- apply(X, 2, function(x) sqrt(n/sum(x^2))*x)

          if (is.null(lambda)) {
            user.lambda<-NULL    #make user.lambda NULL in a call to glmnet
          } else {
            nlambda <- length(lambda)   #override this parameter
            user.lambda <- lambda
          }

          mL <- glmnet::glmnet(Xg, y, alpha = 1, intercept = interc, family = "binomial", nlambda = nlambda, lambda = user.lambda)
          RL <- mL$lambda
          dfy <- apply(mL$beta, 2, function(x) sum(x!=0))
          kt <- 1:length(RL)
          ngp <- which(dfy >= n/2)
          if (length(ngp) > 0){
             RL <- RL[-ngp]
             kt <- kt[-ngp]
             dfy <- dfy[-ngp]
          }
          kk <- which(dfy == 0)
          if(length(kk) > 0){
                        RL <- RL[-kk]
                        kt <- kt[-kk]
          }
          bb <- as.matrix(abs(mL$beta[, kt]))
          SS <- ifelse(bb > 0, 1, 0)
          ii <- duplicated(t(SS))
          bb = bb[,ii == FALSE, drop = FALSE]
          B <- apply(bb, 2, function(x) stats::quantile(x[x!=0], seq(0, 1, length = (o + 1))[-(o + 1)]))
          S <- sapply(1:o, function(j){
            out <- sapply(1:ncol(bb), function(i) ifelse(bb[,i] >= B[j,i], 1, 0))
          })
          SS <- matrix(S, p, sum(ii == FALSE)*o)
          SS <- t(unique(t(SS)))
          mm <- lapply(1:ncol(SS), function(i) SOSnet4glm_help(SS[,i], mL, X, y, lam = lam, interc = interc))
          maxl <- max(sapply(1:length(mm), function(i) length(mm[[i]]$loglik[1,])))
          loglik <- sapply(1:length(mm), function(i) c( unlist(mm[[i]]$loglikbe[1,]), rep(-Inf, maxl - length(unlist(mm[[i]]$loglikbe[1,])))))
          iid <- apply(loglik, 1, which.max)
          if (interc == TRUE){
             maxi <- min(p + 1, maxp)
          } else{
             maxi <- min(p, maxp)
          }
          if (length(iid) > maxi){
             idx <- maxi:1
          } else{
            idx <- length(iid):1
          }
          be = sapply(idx, function(i) {
             return(unlist(mm[[iid[i]]]$loglikbe[2, i]))
          })
          loglik = sapply(idx, function(i) {
             return(loglik = unlist(mm[[iid[i]]]$loglikbe[1, i]))
          })
          fit <- list(beta = be, df = length(idx):1, loglik = loglik, n = n, levels.listed = c(), lambda=mL$lambda, arguments = list(family = "binomial", o = o, nlambda = nlambda, lam = lam, interc = interc, maxp = maxp, lambda = lambda),  interc = interc)
          class(fit) = "DMR"
          return(fit)
}
