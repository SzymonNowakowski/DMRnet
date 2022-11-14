SOSnet4glm <- function(X, y, o, nlambda, lam, interc, maxp, lambda){

          y <- prelasso_binomial(y)

          out <- prelasso_cont_columns(X, y)
          n <-             out$n
          nn <-            out$nn
          p <-             out$p.x + interc     #  in SOSnet p is p.x but maybe with intercept added.
          p.x <-           out$p.x              #  p.x is always without intercept.
          X <-             out$X

          Xg <-            apply(X, 2, function(x) sqrt(n/sum(x^2))*x)

          if(sum(nn != "numeric") > 0){
              stop("Error: wrong data type, columns should be one of types: integer, numeric")
          }

          if (is.null(lambda)) {
            user.lambda<-NULL    #make user.lambda NULL in a call to glmnet
          } else {
            nlambda <- length(lambda)   #override this parameter
            user.lambda <- lambda
          }
###############LASSO#####################
          mL <- glmnet::glmnet(Xg, y, alpha = 1, intercept = interc, family = "binomial", nlambda = nlambda, lambda = user.lambda)
########################################
          bb <- postlasso_common(mL$lambda, n/2, glmnet::coef.glmnet(mL))
          #the calculations were done for mL$beta in v. prior to 0.3.2.9002
          #now, instead of mL$beta (no intercept) I pass coef(mL) which include Intercept. It helps when checks on dfy variable are performed inside
          #also, I pass n/2 instead of n to make it eliminate models with too many variables as in original AP's code

          SS <- postlasso_O_step_preparation(p, p.x, n/2, o, bb[-1, ,drop=FALSE], interc=interc)  # instead of fac we pass bb but without the intercept
                                                                                     # and n/2

          mm <- lapply(1:ncol(SS), function(i) SOSnet4glm_help(SS[,i], mL, X, y, lam = lam, interc = interc))
          maxl <- max(sapply(1:length(mm), function(i) length(mm[[i]]$loglik[1,])))
          loglik <- sapply(1:length(mm), function(i) c( unlist(mm[[i]]$loglikbe[1,]), rep(-Inf, maxl - length(unlist(mm[[i]]$loglikbe[1,])))))

          if (maxl == 1)
            loglik <- t(as.matrix(loglik))      #making loglik a horizontal one-row matrix

          iid <- apply(loglik, 1, which.max)

          maxi <- min(p, maxp)
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
