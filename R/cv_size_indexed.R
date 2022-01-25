cv_size_indexed <- function(X, y, nfolds, model_function, ...) {

        family = list(...)$family

        if (family == "gaussian"){
                n <- length(y)
                real_n <- 0 #recount  of test instances
                foldid <- sample(rep(1:nfolds,length.out=n))   #PP replaces cvfolds by a simpler sample(rep()) function
                error <- list()

                model.full <- model_function(X, y, ...)
                lambda.full<- model.full$lambda

                for (fold in 1:nfolds){
                        Xte <- X[foldid == fold, ,drop = FALSE]
                        yte <- y[foldid == fold]
                        Xtr <- X[foldid != fold, ,drop = FALSE]
                        ytr <- y[foldid != fold]

                        compute_model <- cv_compute_model(model_function, Xtr, ytr, Xte, yte, real_n, lambda.full = lambda.full, ...)   #three letter abbreviations (lambda.full vs lam) make this function call confused, so explicit passing of named parameter i.e. lambda.full=lambda.full is required
                        model<-compute_model$model
                        Xtr<-compute_model$Xtr
                        ytr<-compute_model$ytr
                        Xte<-compute_model$Xte
                        yte<-compute_model$yte
                        real_n<-compute_model$real_n

                        pred <- predict.DMR(model, newx = as.data.frame(Xte))
                        error[[fold]] <- apply(pred, 2, function(z) sum((z - yte)^2))
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

                        foldid1 <- sample(rep(1:nfolds,length.out=n1))  #PP replaces cvfolds by a simpler sample(rep()) function
                        foldid2 <- sample(rep(1:nfolds,length.out=n2))  #PP replaces cvfolds by a simpler sample(rep()) function
                        foldid <- c()
                        foldid[which(y == levels(factor(y))[1])] = foldid1
                        foldid[which(y == levels(factor(y))[2])] = foldid2
                        error <- list()

                        model.full <- model_function(X, y, ...)
                        lambda.full<- model.full$lambda

                        for (fold in 1:nfolds){
                                Xte <- X[foldid == fold, , drop = FALSE]
                                yte <- y[foldid == fold]
                                Xtr <- X[foldid != fold, , drop = FALSE]
                                ytr <- y[foldid != fold]

                                compute_model <- cv_compute_model(model_function, Xtr, ytr, Xte, yte, real_n, lambda.full = lambda.full, ...)   #three letter abbreviations (lambda.full vs lam) make this function call confused, so explicit passing of named parameter i.e. lambda.full=lambda.full is required
                                model<-compute_model$model
                                Xtr<-compute_model$Xtr
                                ytr<-compute_model$ytr
                                Xte<-compute_model$Xte
                                yte<-compute_model$yte
                                real_n<-compute_model$real_n

                                pred <- predict.DMR(model, newx = as.data.frame(Xte), type = "class")
                                error[[fold]] <- apply(pred, 2, function(z) sum(z != yte))
                        }

                }
                else{
                        stop("Error: wrong family, should be one of: gaussian, binomial")
                }
        }
        foldmin <- min(sapply(error, length))
        error <- sapply(1:length(error), function(i) error[[i]][(length(error[[i]]) - foldmin + 1) : length(error[[i]])])
        error <- rowSums(error)/real_n

        kt <- which(error == min(na.omit(error)))
        df.min <- model$df[kt[length(kt)]]

        df.1se = NULL
        kt <- which(error <= min(na.omit(error)) + sd(na.omit(error)))
        df.1se <- model$df[kt[length(kt)]]

        out <- list(df.min = df.min, df.1se = df.1se, dmr.fit = model.full, cvm = error, foldid = foldid)
        return(out)
}
