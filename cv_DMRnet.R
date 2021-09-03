
source("cv_helper.R")

cv_DMRnet <- function(X, y, family = "gaussian", clust.method = 'complete', o = 5, agressive=FALSE, nlambda = 20, lam = 10^(-7), interc = TRUE, nfolds = 10, maxp = ifelse(family == "gaussian", ceiling(length(y)/2), ceiling(length(y)/4))){
        X <- data.frame(X, check.names = TRUE, stringsAsFactors = TRUE)
        if (family == "gaussian"){
                n <- length(y)
                real_n <- 0
                foldid <- sample(rep(1:nfolds,length.out=n))  #same as cvfolds, but written with available functions
                error <- list()
                for (fold in 1:nfolds){
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

                        dmr <- DMRnet(Xtr, ytr, family = "gaussian", clust.method = clust.method, o = o, nlambda = nlambda, interc = interc, maxp = ceiling(maxp))
                        pred <- predict.DMR(dmr, newx = as.data.frame(Xte))
                        error[[fold]] <- apply(pred, 2, function(z) sum((z - yte)^2))
                }
                foldmin <- min(sapply(error, length))
                error <- sapply(1:length(error), function(i) error[[i]][(length(error[[i]]) - foldmin + 1) : length(error[[i]])])
                error <- rowSums(error)/real_n
                dmr.fit <- DMRnet(X, y, family = "gaussian", clust.method = clust.method, o = o, nlambda = nlambda, interc = interc, maxp = ceiling(maxp))
                kt <- which(error == min(error))
                df.min <- dmr$df[kt[length(kt)]]
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
                        real_n <- 0
                        foldid1 <- sample(rep(1:nfolds,length.out=n1))   #same as cvfolds, but written with available functions
                        foldid2 <- sample(rep(1:nfolds,length.out=n2))
                        foldid <- c()
                        foldid[which(y == levels(factor(y))[1])] = foldid1
                        foldid[which(y == levels(factor(y))[2])] = foldid2
                        error <- list()
                        for (fold in 1:nfolds){
                                Xte <- X[foldid == fold, , drop = FALSE]
                                yte <- y[foldid == fold]
                                Xtr <- X[foldid != fold, , drop = FALSE]
                                ytr <- y[foldid != fold]

                                helper <- cv_helper(Xtr, ytr, Xte, yte, real_n, agressive)
                                Xtr<-helper$Xtr
                                ytr<-helper$ytr
                                Xte<-helper$Xte
                                yte<-helper$yte
                                real_n<-helper$real_n

                                dmr <- DMRnet(Xtr, ytr, family = "binomial", clust.method = clust.method, o = o, nlambda = nlambda, lam = lam, interc = interc, maxp = maxp)
                                pred <- predict.DMR(dmr, newx = as.data.frame(Xte), type = "class")
                                error[[fold]] <- apply(pred, 2, function(z) sum(z != yte))
                        }
                        foldmin <- min(sapply(error, length))
                        error <- sapply(1:length(error), function(i) error[[i]][(length(error[[i]]) - foldmin + 1) : length(error[[i]])])
                        error <- rowSums(error)/(real_n)
                        dmr.fit <- DMRnet(X, y, family = "binomial", clust.method = clust.method, o = o, nlambda = nlambda, lam = lam, interc = interc, maxp = maxp)
                        kt <- which(error == min(error))
                        df.min <- dmr$df[kt[length(kt)]]
                }
                else{
                        stop("Error: wrong family, should be one of gaussian, binomial")
                }
        }
        out <- list(df.min = df.min, dmr.fit = dmr.fit, cvm = error, foldid = foldid)
        class(out) <- "cv.DMR"
        return(out)
}
