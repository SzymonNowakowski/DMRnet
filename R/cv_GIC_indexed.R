cv_GIC_indexed <- function(X, y, nfolds, model_function, ...) {

        family = list(...)$family

        if (family == "gaussian"){
                n <- length(y)
                real_n <- 0 #recount  of test instances
                #PP new code foldid <- cvfolds(n, nfolds)
                foldid <- sample(rep(1:nfolds,length.out=n))   #PP replaces cvfolds by a simpler sample(rep()) function
                #PP new code error <- list()
                err <- list(); rss <- list(); #md <- list()

                model.full <- model_function(X, y, ...)
                lambda.full<- model.full$lambda

                for (fold in 1:nfolds){

                        Xte <- X[foldid == fold, ,drop = FALSE]
                        yte <- y[foldid == fold, drop = FALSE]
                        Xtr <- X[foldid != fold, ,drop = FALSE]
                        ytr <- y[foldid != fold, drop = FALSE]

                        compute_model <- cv_compute_model(model_function, Xtr, ytr, Xte, yte, real_n, lambda.full = lambda.full, ...)   #three letter abbreviations (lambda.full vs lam) make this function call confused, so explicit passing of named parameter i.e. lambda.full=lambda.full is required
                        model<-compute_model$model
                        Xtr<-compute_model$Xtr
                        ytr<-compute_model$ytr
                        Xte<-compute_model$Xte
                        yte<-compute_model$yte
                        real_n<-compute_model$real_n

                        #PP new code
                        rss[[fold]] <- model$rss

                        pred <- predict.DMR(model, newx = as.data.frame(Xte))
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




                #PP new code kt <- which(error == min(error))
                #            df.min <- dmr$df[kt[length(kt)]]
                p1 <- model.full$df[1]
                s2 <- model.full$rss[1]/(n-p1)
                Const <- exp(seq(log(2/50),log(2*50), length=80))
                laGIC <- Const*log(p1)*s2
                RSS <- sapply(1:nfolds, function(i) rss[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
                #MD <- sapply(1:nfolds, function(i)  md[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
                IND <- apply( RSS, 2, function(r) sapply( laGIC, function(la) which.min(r+la*length(r):1) ) )
                errGIC <- apply( IND, 1, function(ind) mean(ERR[cbind(ind,1:nfolds)]) )
                #mdGIC  <- apply( IND, 1, function(ind) mean(MD[cbind(ind,1:10)]) )
                #plot(mdGIC[length(laGIC):1],errGIC[length(laGIC):1]/s2, xlab="MD", ylab="PE", type="o")

                r <- model.full$rss
                kt <- which(errGIC == min(errGIC))
                indGIC <- kt[length(kt)]
                gic.full <- (r+laGIC[indGIC]*length(r):1)/(real_n*s2)
                #plot(gic.full[length(gic.full):1])
                indMod <- which.min(gic.full)
                df.min <- model.full$df[indMod]



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
                        real_n <- 0 #recount  of test instances

                        foldid1 <- sample(rep(1:nfolds,length.out=n1))  #PP replaces cvfolds by a simpler sample(rep()) function
                        foldid2 <- sample(rep(1:nfolds,length.out=n2))  #PP replaces cvfolds by a simpler sample(rep()) function
                        foldid <- c()
                        foldid[which(y == levels(factor(y))[1])] = foldid1
                        foldid[which(y == levels(factor(y))[2])] = foldid2
                        #PP new code error <- list()
                        err <- list(); loglik <- list(); #md <- list()

                        model.full <- model_function(X, y, ...)
                        lambda.full<- model.full$lambda

                        for (fold in 1:nfolds) {

                                Xte <- X[foldid == fold, , drop = FALSE]
                                yte <- y[foldid == fold, drop = FALSE]
                                Xtr <- X[foldid != fold, , drop = FALSE]
                                ytr <- y[foldid != fold, drop = FALSE]

                                compute_model <- cv_compute_model(model_function, Xtr, ytr, Xte, yte, real_n, lambda.full = lambda.full, ...)   #three letter abbreviations (lambda.full vs lam) make this function call confused, so explicit passing of named parameter i.e. lambda.full=lambda.full is required
                                model<-compute_model$model
                                Xtr<-compute_model$Xtr
                                ytr<-compute_model$ytr
                                Xte<-compute_model$Xte
                                yte<-compute_model$yte
                                real_n<-compute_model$real_n

                                #SzN new code based on PP new code
                                loglik[[fold]] <- -2*model$loglik

                                pred <- predict.DMR(model, newx = as.data.frame(Xte), type = "class")
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





                        #SzN new code based on PP's new code kt <- which(error == min(error))
                        #            df.min <- dmr$df[kt[length(kt)]]
                        p1 <- model.full$df[1]
                        Const <- exp(seq(log(2/50),log(2*50), length=80))
                        laGIC <- Const*log(p1)
                        LOGLIK <- sapply(1:nfolds, function(i) loglik[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
                        #MD <- sapply(1:nfolds, function(i)  md[[i]][ (len_err[i] - foldmin + 1) : len_err[i] ] )
                        IND <- apply( LOGLIK, 2, function(ll) sapply( laGIC, function(la) which.min(ll+la*length(ll):1) ) )
                        errGIC <- apply( IND, 1, function(ind) mean(ERR[cbind(ind,1:nfolds)]) )
                        #mdGIC  <- apply( IND, 1, function(ind) mean(MD[cbind(ind,1:10)]) )
                        #plot(mdGIC[length(laGIC):1],errGIC[length(laGIC):1]/s2, xlab="MD", ylab="PE", type="o")

                        ll <- -2*model.full$loglik
                        kt <- which(errGIC == min(errGIC))
                        indGIC <- kt[length(kt)]
                        gic.full <- (ll+laGIC[indGIC]*length(ll):1)/real_n
                        #plot(gic.full[length(gic.full):1])
                        indMod <- which.min(gic.full)
                        df.min <- model.full$df[indMod]

                }
                else{
                        stop("Error: wrong family, should be one of: gaussian, binomial")
                }
        }
        #PP: out <- list(df.min = df.min, dmr.fit = dmr.fit, cvm = error, foldid = foldid)
        out <- list(df.min = df.min, df.1se = NULL, dmr.fit = model.full, cvm = gic.full, foldid = foldid)
        return(out)
}
