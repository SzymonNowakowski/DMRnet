cv_MD_indexed <- function(X, y, nfolds, model_function, ...) {

        family = list(...)$family
        algorithm = list(...)$algorithm
        if (!is.null(algorithm))
          if (algorithm == "var_sel")
            warning("Variable selection algorithm was selected. Resulting model sizes may be different in cross validation folds, especially if data with factors with numerous levels is involved or for small cross validation fold counts. This may cause instability in cross validation with model size indexation. If you encounter such problems with your particular data, choose the default GIC indexed cross validation to accompany variable selection algorithm.")

        if (family == "gaussian"){
                n <- length(y)
                total_n <- 0 #recount  of test instances
                foldid <- sample(rep(1:nfolds,length.out=n))   #PP replaces cvfolds by a simpler sample(rep()) function
                error <- list()
                fold_n <- list()

                model.full <- model_function(X, y, ...)
                lambda.full<- model.full$lambda

                for (fold in 1:nfolds){
                        Xte <- X[foldid == fold, ,drop = FALSE]
                        yte <- y[foldid == fold, drop = FALSE]
                        Xtr <- X[foldid != fold, ,drop = FALSE]
                        ytr <- y[foldid != fold, drop = FALSE]

                        compute_model <- cv_compute_model(model_function, Xtr, ytr, Xte, yte, lambda.full = lambda.full, ...)   #three letter abbreviations (lambda.full vs lam) make this function call confused, so explicit passing of named parameter i.e. lambda.full=lambda.full is required
                        model <- compute_model$model
                        Xtr <- compute_model$Xtr
                        ytr <- compute_model$ytr
                        Xte <- compute_model$Xte
                        yte <- compute_model$yte
                        total_n <- total_n + compute_model$this_fold_n
                        fold_n[[fold]] <- compute_model$this_fold_n

                        pred <- predict.DMR(model, newx = as.data.frame(Xte))
                        error[[fold]] <- apply(pred, 2, function(z) sum((z - yte)^2))
                }

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
                        error <- list()

                        model.full <- model_function(X, y, ...)
                        lambda.full<- model.full$lambda

                        for (fold in 1:nfolds){
                                Xte <- X[foldid == fold, , drop = FALSE]
                                yte <- y[foldid == fold, drop = FALSE]
                                Xtr <- X[foldid != fold, , drop = FALSE]
                                ytr <- y[foldid != fold, drop = FALSE]

                                compute_model <- cv_compute_model(model_function, Xtr, ytr, Xte, yte, lambda.full = lambda.full, ...)   #three letter abbreviations (lambda.full vs lam) make this function call confused, so explicit passing of named parameter i.e. lambda.full=lambda.full is required
                                model <- compute_model$model
                                Xtr <- compute_model$Xtr
                                ytr <- compute_model$ytr
                                Xte <- compute_model$Xte
                                yte <- compute_model$yte
                                total_n <- total_n + compute_model$this_fold_n
                                fold_n[[fold]] <- compute_model$this_fold_n

                                pred <- predict.DMR(model, newx = as.data.frame(Xte), type = "class")
                                error[[fold]] <- apply(pred, 2, function(z) sum(z != yte))

                        }

                }
                else{
                        stop("Error: wrong family, should be one of: gaussian, binomial")
                }
        }

        #each fold may have a different length of model sizes (or lambdas)
        #compute foldmin - the length of model sizes (or lambdas) that is minimal across folds.
        foldmin <- min(c(sapply(error, length), length(model.full$df)))   #taking into consideration the length of a full model, which may be SMALLER than in any of the folds

        error_nfolds_length <- length(error[[nfolds]])  #this value needs to be retained because error will be redefined (shifted!) in the next line
        error <- sapply(1:length(error), function(i) error[[i]][(length(error[[i]]) - foldmin + 1) : length(error[[i]])])

        # error is a matrix with foldmin rows and nfolds column,
        # containing sums of errors within each model size and in each fold

        if (foldmin == 1) {
          error <- t(as.matrix(error))  #making it a horizontal one-row matrix
        }

        error_means_in_folds <- error / unlist(fold_n)
        # error_means_in_folds is a matrix with foldmin rows and nfolds column,
        # containing mean errors within each model size and in each fold

        error_means <- rowSums(error)/total_n   # a sum of errors across all folds divided by the total number of test cases
        # error_means is a vector sized foldmin now and it stores mean classification errors for different lambdas (model sizes)
        # note that it is correctly weighted

        wt <- unlist(fold_n) / total_n
        std_err <- sqrt(apply(error_means_in_folds^2-error_means^2, 1, weighted.mean, wt)) / sqrt(nfolds-1)
        # std_err is a vector (of length foldmin) and it stores standard errors for different model sizes
        # calculated similarly as in glmnet for group=TRUE, i.e. doubly weighted (once for mean calculation and then for sd calculation)

        kt <- which(error_means == min(stats::na.omit(error_means)))   #kt stores indexes in error_means equal to a minimum error.
                #if there is more than one such index, the LAST one is the one returned, because LAST means a smaller model.
        df.min <- model$df[error_nfolds_length - foldmin + kt[length(kt)]]   #model is a model computed in the last cross validation fold (==nfolds)
              #so in case there are differences in model lengths in cv folds, the model size in that particular model needs to be shifted


        kt <- which(error_means <= min(stats::na.omit(error_means)) + stats::na.omit(std_err[error_means!=Inf & error_means!=-Inf]))

        if (length(kt) == 0) {
          df.1se <- NULL
        } else {
          df.1se <- model$df[error_nfolds_length - foldmin + kt[length(kt)]]
        }

        out <- list(df.min = df.min, df.1se = df.1se, dmr.fit = model.full, cvm = error_means, cvse = std_err, foldid = foldid)
        return(out)
}
