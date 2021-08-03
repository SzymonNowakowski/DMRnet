

library(randomForest)
library(glmnet)
library(stats)  #glm
library(CatReg)
library(DMRnet)
library(digest)

#library(devtools)
#load_all()

set.seed(strtoi(substr(digest("adult", "md5", serialize = FALSE),1,7),16))



cv_helper<-function(Xtr, ytr, Xte, yte, real_n) {
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


  ####SzN remove from train and test columns causing data singularity
  #preparation to detect singularity
  Xtr.make <- makeX(Xtr)
  prev_pos <- 0
  for (i in 1:ncol(Xtr))
    if (i %in% faki) {  #removing columns from the last level, it is linearly dependant
      # cat(i, prev_pos, length(levels(insurance.train.10percent.x[,i])), "\n")
      Xtr.make <- Xtr.make[,-(prev_pos+length(levels(Xtr[,i])))]
      prev_pos <- prev_pos+length(levels(Xtr[,i])) - 1
    } else prev_pos<-prev_pos+1
  QR<- qr(Xtr.make)

  if (QR$rank < ncol(Xtr.make)) {  #singular
    reverse_lookup<-rep(0, ncol(Xtr.make))
    pos<-1
    for (i in 1:ncol(Xtr))
      if (i %in% faki) {
        reverse_lookup[pos:(pos+length(levels(Xtr[,i]))-1)]<-i  #there are levels-1 columns corresponding to each original column
        pos<-pos+length(levels(Xtr[,i]))-1
      } else {
        reverse_lookup[pos:(pos+1)]<-i
        pos<-pos+1
      }
    #removal of columns for pivot positions larger than rank
    remove_us<-reverse_lookup[QR$pivot[(QR$rank+1):length(QR$pivot)]]
    Xtr <- Xtr[,-remove_us]
    Xte <- Xte[,-remove_us]
    cat("removed", length(unique(remove_us)), "columns\n")
  }
  return (list(Xtr=Xtr, ytr=ytr, Xte=Xte, yte=yte, real_n=real_n))
}


cv_DMRnet <- function(X, y, family = "gaussian", clust.method = 'complete', o = 5, nlambda = 20, lam = 10^(-7), interc = TRUE, nfolds = 10, maxp = ifelse(family == "gaussian", ceiling(length(y)/2), ceiling(length(y)/4))){


       X <- data.frame(X, check.names = TRUE, stringsAsFactors = TRUE)
       if (family == "gaussian"){
          n <- length(y)
          real_n <- 0 #recount  of test instances
          #PP new code foldid <- cvfolds(n, nfolds)
          foldid <- sample(rep(1:nfolds,length.out=n))   #PP replaces nfolds by a simpler sample(rep()) function
          #PP new code error <- list()
          err <- list(); rss <- list(); #md <- list()

          for (fold in 1:nfolds){
              cat("gaussian fold:", fold, "\n")
              Xte <- X[foldid == fold, ,drop = FALSE]
              yte <- y[foldid == fold]
              Xtr <- X[foldid != fold, ,drop = FALSE]
              ytr <- y[foldid != fold]

              helper<- cv_helper(Xtr, ytr, Xte, yte, real_n)
              Xtr<-helper$Xtr
              ytr<-helper$ytr
              Xte<-helper$Xte
              yte<-helper$yte
              real_n<-helper$real_n

              dmr <- DMRnet(Xtr, ytr, family = "gaussian", clust.method = clust.method, o = o, nlambda = nlambda, interc = interc, maxp = ceiling(maxp))
              #PP new code
              rss[[fold]] <- dmr$rss

              pred <- predict(dmr, newx = as.data.frame(Xte))
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
              cat("binomial fold:", fold, "\n")
              Xte <- X[foldid == fold, ,drop = FALSE]
              yte <- y[foldid == fold]
              Xtr <- X[foldid != fold, ,drop = FALSE]
              ytr <- y[foldid != fold]

              helper<- cv_helper(Xtr, ytr, Xte, yte, real_n)
              Xtr<-helper$Xtr
              ytr<-helper$ytr
              Xte<-helper$Xte
              yte<-helper$yte
              real_n<-helper$real_n

              dmr <- DMRnet(Xtr, ytr, family = "binomial", clust.method = clust.method, o = o, nlambda = nlambda, lam = lam, interc = interc, maxp = maxp)
              #SzN new code based on PP new code
              loglik[[fold]] <- -2*dmr$loglik

              pred <- predict(dmr, newx = as.data.frame(Xte), type = "class")
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




adult.train<-read.csv("adult.data", header=FALSE, comment.char="|", stringsAsFactors = TRUE)
adult.test<-read.csv("adult.test", header=FALSE, comment.char="|", stringsAsFactors = TRUE)

colnames <- c('age', #
                     'workclass', #
                     'fnlwgt',
                     'education', #
                     'education_num',
                     'marital_status', #
                     'occupation', #
                     'relationship', #
                     'race', #
                     'sex', #
                     'capital_gain',
                     'capital_loss',
                     'hours_per_week',
                     'native_country', #
                     'income') #
colnames(adult.train)<-colnames
colnames(adult.test)<-colnames
adult.train<-subset(adult.train, adult.train[,14] != " ?")
adult.train<-subset(adult.train, adult.train[,7] != " ?")
adult.train<-subset(adult.train, adult.train[,2] != " ?")
adult.train[,14] <- factor(adult.train[,14])
adult.train[,7] <- factor(adult.train[,7])
adult.train[,2] <- factor(adult.train[,2])
adult.train[,1] <- adult.train[,1] + 0.0  #make a continuous variable out of an integer. Otherwise scope would treat it as a factor
adult.train[,13] <- adult.train[,13] + 0.0

adult.test<-subset(adult.test, adult.test[,14] != " ?")
adult.test<-subset(adult.test, adult.test[,7] != " ?")
adult.test<-subset(adult.test, adult.test[,2] != " ?")
adult.test[,14] <- factor(adult.test[,14])
adult.test[,7] <- factor(adult.test[,7])
adult.test[,2] <- factor(adult.test[,2])
adult.test[,1] <- adult.test[,1] + 0.0  #make a continuous variable out of an integer. Otherwise scope would treat it as a factor
adult.test[,13] <- adult.test[,13] + 0.0

#consiliation of different level names in train and test sets (they end with '.' in test set)
levels(adult.test[,15])[1]<-0
levels(adult.test[,15])[2]<-1
levels(adult.train[,15])[1]<-0
levels(adult.train[,15])[2]<-1


adult.all<-rbind(adult.train, adult.test)
####HURRAY. In total 45222 observations (train+test) as in Stokell's paper

cat("data loaded\n")

errors<-list()
effective_lengths<-list()
sizes<-list()
computation_times<-list()

gamma<-100

#1 PERCENT TRAIN / 99 PERCENT TEST SPLIT
runs<-200
for (model_choice in c( "cv.DMRnet", "gic.DMRnet", "RF", "lr", "cv.glmnet", "scope", "scope")) {
	gamma <- 350 - gamma    #it alternates between 250 and 100
	times<-dfmin<-misclassification_error<-lengths<-rep(0,runs)
	run<-1

	while (run<=runs) {
	  cat("generating train/test sets\n")
	  sample.1percent <- sample(1:nrow(adult.all), 0.01*nrow(adult.all))
	  adult.train.1percent.x <- adult.all[sample.1percent,c(1,2,4,6:10,13:14)] #I exclude education_num and fnlwgt and capital_gain & capital_loss
	  adult.train.1percent.y <- adult.all[sample.1percent,15]

	  adult.test.1percent.x <- adult.all[-sample.1percent,c(1,2,4,6:10,13:14)]
	  adult.test.1percent.y <- adult.all[-sample.1percent,15]


	  #####RECOMPUTATION OF RELEVANT FACTORS in train set, to remove levels with no representative data (empty factors). Needed for random forest and glmnet
	  ###and for DMRnet - old package
	  ###but nor for DMRnet - new package
	  for (i in c(2:8,10))
	    adult.train.1percent.x[,i] <- factor(adult.train.1percent.x[,i])



	  #remove data from test set with factors not present in train subsample as this causes predict() to fail
	  for (i in c(2:8,10)) {
	    train.levels <- levels(adult.train.1percent.x[,i])
	    adult.test.1percent.y<-adult.test.1percent.y[which(adult.test.1percent.x[,i] %in% train.levels)]
	    adult.test.1percent.x<-adult.test.1percent.x[which(adult.test.1percent.x[,i] %in% train.levels),]
	  }
	  for (i in c(2:8,10))
	    adult.test.1percent.x[,i] <- factor(adult.test.1percent.x[,i])   #recalculate factors now for new test


	  start.time <- Sys.time()
	  cat("Started: ", start.time,"\n")

	  if (model_choice=="gic.DMRnet") {
	    cat("DMRnet with GIC only\n")
	    model.1percent <- tryCatch(DMRnet(adult.train.1percent.x, adult.train.1percent.y, nlambda=100, family="binomial"),
	                               error=function(cond) {
	                                 message("Numerical instability in DMRnet detected. Will skip this 1-percent set. Original error:")
	                                 message(cond)
	                                 return(list("red_light"))
	                               })

	    if (model.1percent[[1]] == "red_light") {
	      next
	    }

	    #plot(model.1percent)
	    cat("GIC\n")
	    gic <- gic.DMR(model.1percent, c = 2)
	    #plot(gic)
	  } else  if (model_choice=="cv.DMRnet") {
	      cat("DMRnet with cv\n")
	      model.1percent <- tryCatch(cv_DMRnet(adult.train.1percent.x, adult.train.1percent.y, nlambda=100, family="binomial", nfolds=5),
	                                error=function(cond) {
	                                  message("Numerical instability in cv.DMRnet detected. Will skip this 1-percent set. Original error:")
	                                  message(cond)
	                                  return(list("red_light"))
	                                })

	    if (model.1percent[[1]] == "red_light") {
	      next
	    }

	    #plot(model.1percent)
	    #gic <- gic.DMR(model.1percent, c = 2)
	    #plot(gic)
	  } else if (model_choice=="scope") {
	    cat("Scope, no cv, gamma=", gamma,"\n")
	    model.1percent <- tryCatch(scope.logistic(adult.train.1percent.x, as.numeric(levels(adult.train.1percent.y))[adult.train.1percent.y], gamma=gamma),
	                               error=function(cond) {
	                                 message("Numerical instability in SCOPE detected. Will skip this 1-percent set. Original error:")
	                                 message(cond)
	                                 return(list("red_light"))
	                               })

	    if (model.1percent[[1]] == "red_light") {
	      next
	    }

	  } else if (model_choice=="RF") {
	    cat("random forest. no cv\n")
	    model.1percent <- randomForest(adult.train.1percent.x, y=adult.train.1percent.y)
	  } else if (model_choice=="lr") {
	    cat("Linear Regression no cv\n")
	    model.1percent <- glm(adult.train.1percent.y~., data = adult.train.1percent.x, family="binomial")
	  } else if (model_choice=="cv.glmnet") {
	    cat("glmnet with cv\n")
	    model.1percent<-cv.glmnet(makeX(adult.train.1percent.x), adult.train.1percent.y, family="binomial", nfolds=5)
	  } else
	    stop("Uknown method")




	  if (model_choice=="gic.DMRnet") {
	    cat("DMRnet pred\n")
	    prediction<- tryCatch(predict(model.1percent, newx=adult.test.1percent.x, df = gic$df.min, type="class"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (DMRnet) detected. Will skip this 1-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })

	    if (length(prediction)==2) {
	      next
	    }
	  } else  if (model_choice=="cv.DMRnet") {
	    cat("DMRnet pred\n")
	    prediction<- tryCatch(predict(model.1percent, newx=adult.test.1percent.x, type="class"),#df = gic$df.min, type="class"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (DMRnet) detected. Will skip this 1-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })

	    if (length(prediction)==2) {
	      next
	    }
	  } else if (model_choice=="scope") {
	    cat("scope pred\n")
	    prediction<- ifelse(predict(model.1percent, adult.test.1percent.x) >0.5,1,0)
	  } else if (model_choice=="RF") {
	    cat("Random Forest pred\n")
	    prediction<- tryCatch(predict(model.1percent, adult.test.1percent.x, type="class"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (RF) detected. Will skip this 1-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })

	    if (length(prediction)==2) {
	      next
	    }
	  } else if (model_choice=="lr") {
	    cat("Linear Regression pred\n")
	    prediction<- ifelse(predict(model.1percent, adult.test.1percent.x) >0,1,0)
	  } else if (model_choice=="cv.glmnet") {
	    cat("glmnet pred\n")
	    prediction<- tryCatch(predict(model.1percent, newx=makeX(adult.test.1percent.x), type="class"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (cv.glmnet) detected. Will skip this 1-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })

	    if (length(prediction)==2) {
	      next
	    }
	  } else
	    stop("Uknown method")

	  end.time <- Sys.time()
	  times[run] <- as.numeric(end.time)-as.numeric(start.time)
	  cat("Ended: ", end.time,"elapsed: ", times[run],"\n")


	  lengths[run]<-length(prediction[!is.na(prediction)])

	  prediction[is.na(prediction)] <- 0
	  misclassification_error[run]<-mean(prediction[!is.na(prediction)] != adult.test.1percent.y[!is.na(prediction)])

	  if (model_choice == "gic.DMRnet")
	    dfmin[run]<-gic$df.min
	  if (model_choice == "cv.DMRnet" )
	    dfmin[run]<-model.1percent$df.min
	  if (model_choice == "scope")
	    dfmin[run]<-length(unique(c(sapply(sapply(model.1percent$beta.best[[2]], as.factor), levels), sapply(sapply(model.1percent$beta.best[[1]], as.factor), levels),recursive=TRUE)))-1 + #-1 is for "0" level
	                -sum(sapply(sapply(model.1percent$beta.best[[2]], as.factor), levels)!="0")   #and we subtract the number of factors = number of constraints from eq. (8) in Stokell et al.

	  cat(run, "median = ", median(misclassification_error[misclassification_error>0]), "\n")
	  cat(run, "df.min = ", mean(dfmin[misclassification_error>0]), "\n")
	  cat(run, "lengths = ", mean(lengths[misclassification_error>0]), "\n")

	  run<-run+1
	}

	cat("overall median = ", median(misclassification_error[misclassification_error!=0]), "\n")


	model_name<-model_choice
	if (model_choice == "scope")
	  model_name<-paste(model_name, gamma, sep="-")


	computation_times[[model_name]]<-times
	effective_lengths[[model_name]]<-lengths
	if (length(dfmin[dfmin>0])>0)
		sizes[[model_name]]<-dfmin
	errors[[model_name]]<-misclassification_error

}

write.csv(errors, "adult_errors.csv")
write.csv(effective_lengths, "adult_effective_lengths.csv")
write.csv(sizes, "adult_model_sizes.csv")
write.csv(computation_times, "adult_computation_times.csv")


pdf("adult_computation_times.pdf",width=12,height=5)
boxplot(computation_times)
dev.off

pdf("adult_errors.pdf",width=12,height=5)
boxplot(errors, ylim=c(0.16, 0.26))
dev.off

pdf("adult_model_sizes.pdf",width=9,height=5)
boxplot(sizes)
dev.off

pdf("adult_effective_lengths.pdf",width=12,height=5)
boxplot(effective_lengths)
dev.off


