

library(randomForest)
library(glmnet)
library(stats)  #glm
library(CatReg)
library(DMRnet)
library(digest)
#
# library(devtools)
# load_all()

set.seed(strtoi(substr(digest("insurance", "md5", serialize = FALSE),1,7),16))


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
        reverse_lookup[pos:(pos+length(levels(Xtr[,i]))-2)]<-i  #there are levels-1 columns corresponding to each original column
        pos<-pos+length(levels(Xtr[,i]))-1
      } else {
        reverse_lookup[pos]<-i
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


insurance.all<-read.csv("insurance_train", header=TRUE, comment.char="|", stringsAsFactors = TRUE)
insurance.all<-insurance.all[,apply(apply(insurance.all,2,is.na), 2, sum)==0]  #removing columns with NA
insurance.all.x<-insurance.all[,2:(ncol(insurance.all)-1)]  #without ID and the response columns
cont_columns = c(4, 8, 9, 10, 11)
number_of_levels<-rep(0, ncol(insurance.all.x))
for (i in 1:ncol(insurance.all.x))
  if (!(i %in% cont_columns)) {
    insurance.all.x[,i] <- factor(insurance.all.x[,i])  #int->factors
    number_of_levels[i] <- length(levels(insurance.all.x[,i]))
  }

cat("data loaded\n")

errors<-list()
effective_lengths<-list()
sizes<-list()
computation_times<-list()

gamma<-8

#1 PERCENT TRAIN / 99 PERCENT TEST SPLIT
runs<-5
for (model_choice in c(  "cv.DMRnet", "gic.DMRnet", "lr",  "scope", "scope")) {
	gamma <- 40 - gamma    #it alternates between 32 and 8
	times<-dfmin<-MSPE<-lengths<-rep(0,runs)
	run<-1

	while (run<=runs) {
	  cat("generating response vector\n")

	  response_factor_columns <- sample(1:length(number_of_levels), 10, prob = number_of_levels/sum(number_of_levels))
	              #in this sample we include continous columns but with 0 probability
	  K<-number_of_levels[response_factor_columns]
	  s<-floor(2+0.5*log(K))
	  theta_coefficients<-list()
	  for (i in 1:10) {
	    theta_coefficients[[i]]<-c(0) #min(s) = 2
	    while (length(unique(theta_coefficients[[i]])) < s[i])  #we sample until the condition "yielding sj true levels" is satisfied
	      theta_coefficients[[i]] <- sample(1:s[i], size=K[i], replace = TRUE)
	  }

	  continous_coefficients <- rnorm(5, 0, 1)

	  insurance.all.y <- rep(0, nrow(insurance.all.x))
	  for (i in 1:10)
	    insurance.all.y <- insurance.all.y + theta_coefficients[[i]][as.integer(insurance.all.x[,response_factor_columns[i]])]
	  for (i in 1:5)
	    insurance.all.y <- insurance.all.y + continous_coefficients[i]*insurance.all.x[,cont_columns[i]]

	  m<-mean(insurance.all.y)
	  std<-sd(insurance.all.y)

	  insurance.all.y <- (insurance.all.y - m)/std+m + rnorm(1, 0, 1) #response was then scaled to have unit variance, after which standard normal noise was added.

	  cat("generating train/test sets\n")

	  sample.10percent <- sample(1:nrow(insurance.all.x), 0.1*nrow(insurance.all.x))
    insurance.train.10percent.x <- insurance.all.x[sample.10percent,]
	  insurance.train.10percent.y <- insurance.all.y[sample.10percent]
	  insurance.test.10percent.x <- insurance.all.x[-sample.10percent,]
	  insurance.test.10percent.y <- insurance.all.y[-sample.10percent]



	  #####RECOMPUTATION OF RELEVANT FACTORS in train set, to remove levels with no representative data (empty factors). Needed for random forest and glmnet
	  ###and for DMRnet - old package
	  ###but nor for DMRnet - new package


	  for (i in 1:ncol(insurance.train.10percent.x))
	    if (!(i %in% cont_columns))
	      insurance.train.10percent.x[,i] <- factor(insurance.train.10percent.x[,i])


	  #remove data from test set with factors not present in train subsample as this causes predict() to fail
	  for (i in 1:ncol(insurance.train.10percent.x))
	    if (!(i %in% cont_columns)) {
	      train.levels <- levels(insurance.train.10percent.x[,i])
	      insurance.test.10percent.y<-insurance.test.10percent.y[which(insurance.test.10percent.x[,i] %in% train.levels)]
	      insurance.test.10percent.x<-insurance.test.10percent.x[which(insurance.test.10percent.x[,i] %in% train.levels),]
	    }
	  for (i in 1:ncol(insurance.train.10percent.x))
	    if (!(i %in% cont_columns)) {
	      train.levels <- levels(insurance.train.10percent.x[,i])
	      insurance.test.10percent.x[,i] <- factor(insurance.test.10percent.x[,i], levels=train.levels)   #recalculate factors now for new test
	    }

	  #handling the singular case
	  insurance.make <- makeX(insurance.train.10percent.x)
	  prev_pos <- 0
	  for (i in 1:ncol(insurance.train.10percent.x))
	    if (!(i %in% cont_columns)) {  #removing columns from the last level, it is linearly dependant
	     # cat(i, prev_pos, length(levels(insurance.train.10percent.x[,i])), "\n")
	      insurance.make<-insurance.make[,-(prev_pos+length(levels(insurance.train.10percent.x[,i])))]
	      prev_pos <- prev_pos+length(levels(insurance.train.10percent.x[,i])) - 1
	    } else prev_pos<-prev_pos+1
	  QR<- qr(insurance.make)

	  if (QR$rank < ncol(insurance.make)) {  #singular
	    reverse_lookup<-rep(0, ncol(insurance.make))
	    pos<-1
	    for (i in 1:ncol(insurance.train.10percent.x))
	      if (!(i %in% cont_columns)) {
	        reverse_lookup[pos:(pos+length(levels(insurance.train.10percent.x[,i]))-2)]<-i  #there are levels-1 columns corresponding to each original column
	        pos<-pos+length(levels(insurance.train.10percent.x[,i]))-1
	      } else {
	        reverse_lookup[pos]<-i
	        pos<-pos+1
	      }
	    #removal of columns for pivot positions larger than rank
	    remove_us<-reverse_lookup[QR$pivot[(QR$rank+1):length(QR$pivot)]]
	    insurance.train.10percent.x <- insurance.train.10percent.x[,-remove_us]
	    insurance.test.10percent.x <- insurance.test.10percent.x[,-remove_us]
	    cat("removed", length(unique(remove_us)), "columns\n")
	  }

	  cat("consolidated factors and columns\n")

	  start.time <- Sys.time()
	  cat("Started: ", start.time,"\n")

	  if (model_choice=="gic.DMRnet") {
	    cat("DMRnet with GIC only\n")
	    model.10percent <- tryCatch(DMRnet(insurance.train.10percent.x, insurance.train.10percent.y, nlambda=100, family="gaussian"),
	                               error=function(cond) {
	                                 message("Numerical instability in DMRnet detected. Will skip this 10-percent set. Original error:")
	                                 message(cond)
	                                 return(list("red_light"))
	                               })

	    if (model.10percent[[1]] == "red_light") {
	      next
	    }

	    #plot(model.10percent)
	    cat("GIC\n")
	    gic <- gic.DMR(model.10percent, c = 2)
	    #plot(gic)
	  } else  if (model_choice=="cv.DMRnet") {
	      cat("DMRnet with cv\n")
	      model.10percent <- tryCatch(cv_DMRnet(insurance.train.10percent.x, insurance.train.10percent.y, nlambda=100, family="gaussian", nfolds=5),
	                                error=function(cond) {
	                                  message("Numerical instability in cv.DMRnet detected. Will skip this 10-percent set. Original error:")
	                                  message(cond)
	                                  return(list("red_light"))
	                                })

	    if (model.10percent[[1]] == "red_light") {
	      next
	    }

	    #plot(model.10percent)
	    #gic <- gic.DMR(model.10percent, c = 2)
	    #plot(gic)
	  } else if (model_choice=="scope") {
	    cat("Scope, no cv, gamma=", gamma,"\n")
	    model.10percent <- tryCatch(scope(insurance.train.10percent.x, insurance.train.10percent.y, gamma=gamma),
	                               error=function(cond) {
	                                 message("Numerical instability in SCOPE detected. Will skip this 10-percent set. Original error:")
	                                 message(cond)
	                                 return(list("red_light"))
	                               })

	    if (model.10percent[[1]] == "red_light") {
	      next
	    }

	  } else if (model_choice=="RF") {
	    cat("random forest. no cv\n")
	    model.10percent <- randomForest(insurance.train.10percent.x, y=insurance.train.10percent.y)
	  } else if (model_choice=="lr") {
	    cat("Linear Regression no cv\n")
	    model.10percent <- glm(insurance.train.10percent.y~., data = insurance.train.10percent.x, family="gaussian")
	  } else if (model_choice=="cv.glmnet") {
	    cat("glmnet with cv\n")
	    model.10percent<-cv.glmnet(makeX(insurance.train.10percent.x), insurance.train.10percent.y, family="gaussian", nfolds=5)
	  } else
	    stop("Uknown method")




	  if (model_choice=="gic.DMRnet") {
	    cat("DMRnet pred\n")
	    prediction<- tryCatch(predict(model.10percent, newx=insurance.test.10percent.x, df = gic$df.min, type="response"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (DMRnet) detected. Will skip this 10-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })

	    if (length(prediction)==2) {
	      next
	    }
	  } else  if (model_choice=="cv.DMRnet") {
	    cat("DMRnet pred\n")
	    prediction<- tryCatch(predict(model.10percent, newx=insurance.test.10percent.x, type="response"),#df = gic$df.min, type="class"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (DMRnet) detected. Will skip this 10-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })

	    if (length(prediction)==2) {
	      next
	    }
	  } else if (model_choice=="scope") {
	    cat("scope pred\n")
	    prediction<- predict(model.10percent, insurance.test.10percent.x)
	  } else if (model_choice=="RF") {
	    cat("Random Forest pred\n")
	    prediction<- tryCatch(predict(model.10percent, insurance.test.10percent.x, type="response"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (RF) detected. Will skip this 10-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })

	    if (length(prediction)==2) {
	      next
	    }
	  } else if (model_choice=="lr") {
	    cat("Linear Regression pred\n")
	    prediction<- predict(model.10percent, insurance.test.10percent.x)
	  } else if (model_choice=="cv.glmnet") {
	    cat("glmnet pred\n")
	    prediction<- tryCatch(predict(model.10percent, newx=makeX(insurance.test.10percent.x), type="response"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (cv.glmnet) detected. Will skip this 10-percent set. Original error:")
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

	  MSPE[run]<-mean((prediction[!is.na(prediction)] - insurance.test.10percent.y[!is.na(prediction)])^2)

	  if (model_choice == "gic.DMRnet")
	    dfmin[run]<-gic$df.min
	  if (model_choice == "cv.DMRnet" )
	    dfmin[run]<-model.10percent$df.min
	  if (model_choice == "scope")
	    dfmin[run]<-length(unique(c(sapply(sapply(model.10percent$beta.best[[2]], as.factor), levels), sapply(sapply(model.10percent$beta.best[[1]], as.factor), levels),recursive=TRUE)))-1 + #-1 is for "0" level
	                -sum(sapply(sapply(model.10percent$beta.best[[2]], as.factor), levels)!="0")   #and we subtract the number of factors = number of constraints from eq. (8) in Stokell et al.

	  cat(run, "median = ", median(MSPE[MSPE>0]), "\n")
	  cat(run, "df.min = ", mean(dfmin[MSPE>0]), "\n")
	  cat(run, "lengths = ", mean(lengths[MSPE>0]), "\n")

	  run<-run+1
	}

	cat("overall median = ", median(MSPE[MSPE!=0]), "\n")


	model_name<-model_choice
	if (model_choice == "scope")
	  model_name<-paste(model_name, gamma, sep="-")


	computation_times[[model_name]]<-times
	effective_lengths[[model_name]]<-lengths
	if (length(dfmin[dfmin>0])>0)
		sizes[[model_name]]<-dfmin
	errors[[model_name]]<-MSPE

}

write.csv(errors, "insurance_errors.csv")
write.csv(effective_lengths, "insurance_effective_lengths.csv")
write.csv(sizes, "insurance_model_sizes.csv")
write.csv(computation_times, "insurance_computation_times.csv")


pdf("insurance_computation_times.pdf",width=12,height=5)
boxplot(computation_times)
dev.off

pdf("insurance_errors.pdf",width=12,height=5)
boxplot(errors, ylim=c(0.16, 0.26))
dev.off

pdf("insurance_model_sizes.pdf",width=9,height=5)
boxplot(sizes)
dev.off

pdf("insurance_effective_lengths.pdf",width=12,height=5)
boxplot(effective_lengths)
dev.off


