SOSnet4lm <- function(X, y, o, nlambda, interc, maxp, lambda){
#SOSnet algorithm pseudocode available from https://arxiv.org/pdf/1907.03025v1.pdf page 17

          out <- prelasso_cont_columns(X, y)
          n <-             out$n
          nn <-            out$nn
          p <-             out$p.x + interc     #  in SOSnet p is p.x but maybe with intercept added.
          p.x <-           out$p.x              #  p.x is always without intercept.

          if(sum(nn != "numeric") > 0){
            stop("Error: wrong data type, columns should be one of types: integer, numeric")
          }

          Xg <- apply(X, 2, function(x) sqrt(n/sum(x^2))*x)

          if (is.null(lambda)) {
            user.lambda<-NULL    #make user.lambda NULL in a call to glmnet
          } else {
            nlambda <- length(lambda)   #override this parameter
            user.lambda <- lambda
          }

###############LASSO#####################
          mL <- glmnet::glmnet(Xg, y, alpha = 1, intercept = interc, family = "gaussian", nlambda = nlambda, lambda = user.lambda)
########################################
          bb <- postlasso_common(mL$lambda, n, coef(mL))
          #the calculations were done for mL$beta in v. prior to 0.3.2.9002
          #now, instead of mL$beta (no intercept) I pass coef(mL) which include Intercept. It helps when checks on dfy variable are performed inside

          SS <- postlasso_O_step_preparation(p, p.x, n, o, bb[-1, ,drop=FALSE], interc=interc)  #instead of fac we pass bb but without the intercept

          mm <- lapply(1:ncol(SS), function(i) SOSnet4lm_help(SS[,i], mL, X, y, interc = interc))

          #at that point we can have duplicated indices too as a result of considering matrices of non-full rank
          duplicated_indices <- duplicated(sapply(1:length(mm), function(i) mm[[i]]$ind))
          mm<-mm[!duplicated_indices]  #removing the duplicates

          maxl <- max(sapply(1:length(mm), function(i) length(mm[[i]]$rss)))
          rss <- sapply(1:length(mm), function(i) c(rep(Inf, maxl - length(mm[[i]]$rss)), mm[[i]]$rss))

          if (maxl == 1)
            rss <- t(as.matrix(rss))      #making rss a horizontal one-row matrix

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
          fit <- list(beta = be, df = length(idx):1, rss = rss[cbind(idx, iid[idx])], n = n, levels.listed = c(), lambda=mL$lambda, arguments = list(family = "gaussian", nlambda = nlambda, interc = interc, maxp = maxp, lambda = lambda), interc = interc)
          class(fit) = "DMR"
          return(fit)
}
