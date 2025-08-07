cv_GIC_indexed <- function(X, y, nfolds, model_function, ...) {

        family = list(...)$family

        if (family == "gaussian"){
                n <- length(y)
                total_n <- 0 #recount  of test instances
                foldid <- sample(rep(1:nfolds,length.out=n))   #PP replaces cvfolds by a simpler sample(rep()) function
                err <- list()
                rss <- list() #md <- list()
                fold_n <-list()

                model.full <- model_function(X, y, ...)
                lambda.full<- model.full$lambda

                for (fold in 1:nfolds){

                        Xte <- X[foldid == fold, ,drop = FALSE]
                        yte <- y[foldid == fold, drop = FALSE]
                        Xtr <- X[foldid != fold, ,drop = FALSE]
                        ytr <- y[foldid != fold, drop = FALSE]

                        compute_model <- cv_compute_model(model_function, Xtr, ytr, Xte, yte, lambda.full = lambda.full, ...)   #three letter abbreviations (lambda.full vs lam) make this function call confused, so explicit passing of named parameter i.e. lambda.full=lambda.full is required
                        model<-compute_model$model
                        Xtr<-compute_model$Xtr
                        ytr<-compute_model$ytr
                        Xte<-compute_model$Xte
                        yte<-compute_model$yte
                        total_n <- total_n + compute_model$this_fold_n
                        fold_n[[fold]] <- compute_model$this_fold_n

                        #PP new code
                        rss[[fold]] <- model$rss

                        pred <- predict.DMR(model, newx = as.data.frame(Xte))
                        #PP new code error[[fold]] <- apply(pred, 2, function(z) sum((z - yte)^2))
                        err[[fold]] <- apply(pred, 2, function(z) mean((z - yte)^2))

                }


                len_err <- sapply(err, length)
                foldmin <- min(len_err)
                ERR <- sapply(1:nfolds, function(i) err[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
                if (foldmin == 1) {
                  ERR<-t(as.matrix((ERR)))  #making it a horizontal one-row matrix
                }
                #err <- rowMeans(ERR); kt <- which(err == min(err)); df.min <- dmr$df[kt[length(kt)]]; plot(err, type="o")

                p1 <- model.full$df[1]
                s2 <- model.full$rss[1]/(n-p1)

                p <- ncol(model.full$beta)
                if (is.null(p))
                  p <- length(model.full$beta)

                RIC_constant <- constants()$RIC_gaussian_constant
                Const <- exp(seq(log(RIC_constant/50),log(RIC_constant*50), length=81))
                laGIC <- Const*log(p)*s2
                RSS <- sapply(1:nfolds, function(i) rss[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
                if (is.null(dim(RSS))) {
                  RSS<-t(as.matrix((RSS)))  #making it a horizontal one-row matrix
                }
                #MD <- sapply(1:nfolds, function(i)  md[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
                IND <- apply( RSS, 2, function(r) sapply( laGIC, function(la) which.min(r+la*length(r):1) ) )
                error_means_in_folds <- apply( IND, 1, function(ind) ERR[cbind(ind,1:nfolds)])
                #here the errors are mean errors for each model size and each fold, sized foldmin x nfolds

                errGIC <- apply(error_means_in_folds, 1, mean)
                #errGIC is the mean of those mean errors


                std_err <- sqrt(rowSums(error_means_in_folds^2-errGIC^2)/((nfolds-1)*nfolds))     # SD in nfolds numbers, further divided by sqrt(n_folds)
                # std_err is a vector (of length foldmin) and it stores standard errors for different model sizes
                # it is calculated similarly as in glmnet for group=TRUE, but not doubly weighted

                #mdGIC  <- apply( IND, 1, function(ind) mean(MD[cbind(ind,1:10)]) )
                #plot(mdGIC[length(laGIC):1],errGIC[length(laGIC):1]/s2, xlab="MD", ylab="PE", type="o")

                r <- model.full$rss
                kt <- which(errGIC == min(errGIC))
                indGIC <- kt[length(kt)]    #TODO: why last?
                gic.full <- (r+laGIC[indGIC]*length(r):1)/(total_n*s2)
                #plot(gic.full[length(gic.full):1])

        } else{
                if (family == "binomial"){
                        if (!inherits(y, "factor")){
                                stop("Error: y should be a factor")
                        }
                        lev <- levels(factor(y))
                        if (length(lev) != 2){
                                stop("Error: factor y should have 2 levels")
                        }
                        n1 <- table(y)[1]
                        n2 <- table(y)[2]
                        total_n <- 0 #recount  of test instances

                        foldid1 <- sample(rep(1:nfolds,length.out=n1))  #PP replaces cvfolds by a simpler sample(rep()) function
                        foldid2 <- sample(rep(1:nfolds,length.out=n2))  #PP replaces cvfolds by a simpler sample(rep()) function
                        foldid <- c()
                        foldid[which(y == levels(factor(y))[1])] = foldid1
                        foldid[which(y == levels(factor(y))[2])] = foldid2
                        #PP new code error <- list()
                        err <- list()
                        loglik <- list()
                        #md <- list()
                        fold_n <-list()

                        model.full <- model_function(X, y, ...)
                        lambda.full<- model.full$lambda

                        for (fold in 1:nfolds) {

                                Xte <- X[foldid == fold, , drop = FALSE]
                                yte <- y[foldid == fold, drop = FALSE]
                                Xtr <- X[foldid != fold, , drop = FALSE]
                                ytr <- y[foldid != fold, drop = FALSE]

                                compute_model <- cv_compute_model(model_function, Xtr, ytr, Xte, yte, lambda.full = lambda.full, ...)   #three letter abbreviations (lambda.full vs lam) make this function call confused, so explicit passing of named parameter i.e. lambda.full=lambda.full is required
                                model<-compute_model$model
                                Xtr<-compute_model$Xtr
                                ytr<-compute_model$ytr
                                Xte<-compute_model$Xte
                                yte<-compute_model$yte
                                total_n <- total_n + compute_model$this_fold_n
                                fold_n[[fold]] <- compute_model$this_fold_n


                                #SzN new code based on PP new code
                                loglik[[fold]] <- -2*model$loglik

                                pred <- predict.DMR(model, newx = as.data.frame(Xte), type = "class")
                                #SzN new code based on PP new code error[[fold]] <- apply(pred, 2, function(z) sum(z != yte))
                                err[[fold]] <- apply(pred, 2, function(z) mean(z != yte))

                        }

                        len_err <- sapply(err, length)
                        foldmin <- min(len_err)
                        ERR <- sapply(1:nfolds, function(i) err[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
                        #err <- rowMeans(ERR); kt <- which(err == min(err)); df.min <- dmr$df[kt[length(kt)]]; plot(err, type="o")


                        if (foldmin == 1) {
                          ERR<-t(as.matrix((ERR)))  #making it a horizontal one-row matrix
                        }

                        p <- ncol(model.full$beta)
                        if (is.null(p))
                          p <- length(model.full$beta)

                        RIC_constant <- constants()$RIC_binomial_constant
                        Const <- exp(seq(log(RIC_constant/50),log(RIC_constant*50), length=81))
                        laGIC <- Const*log(p)
                        LOGLIK <- sapply(1:nfolds, function(i) loglik[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
                        if (is.null(dim(LOGLIK))) {
                          LOGLIK<-t(as.matrix((LOGLIK)))  #making it a horizontal one-row matrix
                        }
                        #MD <- sapply(1:nfolds, function(i)  md[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
                        IND <- apply( LOGLIK, 2, function(ll) sapply( laGIC, function(la) which.min(ll+la*length(ll):1) ) )
                        error_means_in_folds <- apply( IND, 1, function(ind) ERR[cbind(ind,1:nfolds)])
                        #here the errors are mean errors for each model size and each fold, sized foldmin x nfolds

                        errGIC <- apply(error_means_in_folds, 1, mean)
                        #errGIC is the mean of those mean errors


                        std_err <- sqrt(rowSums(error_means_in_folds^2-errGIC^2)/((nfolds-1)*nfolds))     # SD in nfolds numbers, further divided by sqrt(n_folds)
                        # std_err is a vector (of length foldmin) and it stores standard errors for different model sizes
                        # it is calculated similarly as in glmnet for group=TRUE, but not doubly weighted

                        #mdGIC  <- apply( IND, 1, function(ind) mean(MD[cbind(ind,1:10)]) )
                        #plot(mdGIC[length(laGIC):1],errGIC[length(laGIC):1]/s2, xlab="MD", ylab="PE", type="o")

                        ll <- -2*model.full$loglik
                        kt <- which(errGIC == min(errGIC))
                        indGIC <- kt[length(kt)]    #TODO: why last?
                        gic.full <- (ll+laGIC[indGIC]*length(ll):1)/total_n
                        #plot(gic.full[length(gic.full):1])

                }
                else{
                        stop("Error: wrong family, should be one of: gaussian, binomial")
                }
        }
        kt <- which(gic.full == min(stats::na.omit(gic.full))) #kt stores indexes in error equal to a minimum error.
        #if there is more than one such index, the LAST one is the one returned, because LAST means a smaller model.
        indMod <- kt[length(kt)]
        df.min <- model.full$df[indMod]

        kt <- which(gic.full <= min(stats::na.omit(gic.full)) + stats::na.omit(std_err[gic.full!=Inf & gic.full!=-Inf]))
        if (length(kt) == 0) {
          df.1se <- NULL
        } else {
          indMod <- kt[length(kt)]
          df.1se <- model.full$df[indMod]
        }

        out <- list(df.min = df.min, df.1se = df.1se, dmr.fit = model.full, cvm = gic.full, cvse = std_err, foldid = foldid)
        return(out)
}
