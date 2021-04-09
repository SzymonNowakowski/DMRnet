SOSnet4lm <- function(X, y, o = 5, nlambda = 50, interc = TRUE, maxp = ceiling(length(y)/2)){
#SOSnet algorithm pseudocode available from https://arxiv.org/pdf/1907.03025v1.pdf page 17
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
          mL <- glmnet::glmnet(Xg, y, alpha = 1, intercept = interc, nlambda = nlambda, family = "gaussian")
          RL <- mL$lambda
          dfy <- apply(mL$beta, 2, function(x) sum(x!=0))  #equal to s_k from SOSnet pseudocode
          kt <- 1:length(RL)    #indices of lambdas
          ngp <- which(dfy >= n)   #(1) removing predictor sets with more predictors than matrix rows
          if (length(ngp) > 0){
            RL <- RL[-ngp]  #removing them from lambdas
            kt <- kt[-ngp]   #and from lambda indices
            dfy <- dfy[-ngp]
          }
          kk <- which(dfy == 0)    #(2) removing predictor sets with 0 predictors
          if(length(kk) > 0){
            RL <- RL[-kk]     #removing them from lambdas
            kt <- kt[-kk]     #and for lambda indices
          }
          bb <- as.matrix(abs(mL$beta[, kt]))    #Betas corresponding to legal lambdas
          SS <- ifelse(bb > 0, 1, 0)
          ii <- duplicated(t(SS))    #detecting duplicated predictor sets
          bb <- bb[, ii == FALSE, drop = FALSE]    #(3) removing duplicated predictor sets
          B <- apply(bb, 2, function(x) stats::quantile(x[x!=0], seq(0, 1, length = (o + 1))[-(o + 1)]))
                                     #o-based quantiles of non-zero Betas (without the last, 100% quantile)
  ##now, perform the selection of s_kl = floor(s_k*l/o) highest Betas from the k-th step and l-th substep of the SOSnet pseudocode
          #first, note that sum(ii==FALSE) is a number of predictor sets and it may be smaller than nlambda because of (1), (2), (3)
          S <- sapply(1:o, function(j){
            sapply(1:ncol(bb), function(i) ifelse(bb[,i] >= B[j,i], 1, 0))
          })
          SS <- matrix(S, p, sum(ii == FALSE)*o)
  ## a column number (l-1)*sum(ii==FALSE) + k of SS contains a "1" in a given row iff a given highest beta is a part of J i.e. that predictor is a part of s_kl selected predictors
          SS <- t(unique(t(SS)))   #again, removing duplicate rows
          mm <- lapply(1:ncol(SS), function(i) SOSnet4lm_help(i, SS[,i], mL, X, y, interc = interc))

          #at that point we can have duplicated indices too as a result of considering matrices of non-full rank
          duplicated_indices <- duplicated(sapply(1:length(mm), function(i) mm[[i]]$ind))
          mm<-mm[!duplicated_indices]  #removing the duplicates

          maxl <- max(sapply(1:length(mm), function(i) length(mm[[i]]$rss)))
          rss <- sapply(1:length(mm), function(i) c(rep(Inf, maxl - length(mm[[i]]$rss)), mm[[i]]$rss))
          iid <- apply(rss, 1, which.min)

          maxi <- min(p, maxp)
          if (length(iid) > maxi){
             idx <- (length(iid) - maxi):length(iid)
          } else{
            idx <- 1:length(iid)
          }
          if (interc == FALSE){
            be = sapply(idx, function(i) {
              Xs <- X[, mm[[iid[i]]]$ind[1:(length(iid) + 1 - i)], drop = FALSE]
              mnk <- stats::lm.fit(as.matrix(Xs), y)
              out <- rep(0, p)
              out[mm[[iid[i]]]$ind[1:(length(iid) + 1 - i)]] <- mnk$coef
              return(out)
            })
          } else{
             be = sapply(idx[-length(idx)], function(i) {
              Xs <- X[, mm[[iid[i]]]$ind[1:(length(iid) - i)], drop = FALSE]
              mnk <- stats::lm.fit(as.matrix(cbind(1, Xs)), y)
              out <- rep(0, p + 1)
              out[c(1, mm[[iid[i]]]$ind[1:(length(iid) - i)] + 1)] <- mnk$coef
              return(out)
            })
            mnk <- stats::lm.fit(as.matrix(rep(1, n)), y)
            be <- cbind(be, c(mnk$coef, rep(0, ncol(X))))
          }
          fit <- list(beta = be, df = length(idx):1, rss = rss[cbind(idx, iid[idx])], n = n, arguments = list(family = "gaussian", nlambda = nlambda, interc = interc, maxp = maxp), interc = interc)
          class(fit) = "DMR"
          return(fit)
}
