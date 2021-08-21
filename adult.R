

library(randomForest)
library(glmnet)
library(stats)  #glm
library(CatReg)
library(DMRnet)
library(digest)

#library(devtools)
#load_all()

set.seed(strtoi(substr(digest("adult", "md5", serialize = FALSE),1,7),16))

source("glaf_4glm.R")
source("cv_DMRnet.R")

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
runs<-1000
for (model_choice in c( "cv.GLAF", "gic.GLAF", "cv.DMRnet", "gic.DMRnet", "RF", "lr", "cv.glmnet", "scope", "scope")) {
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


	  #removing columns with only one level:
	  singular_factors<-which(sapply(sapply(adult.train.1percent.x, levels), length)==1)   #for continous columns length is 0
	  if (length(singular_factors)>0) {
  	  adult.test.1percent.x <- adult.test.1percent.x[,-singular_factors]
  	  adult.train.1percent.x <- adult.train.1percent.x[,-singular_factors]
  	  cat("removed", length(singular_factors), "columns due to singular factors\n")
	  }

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

	    cat("GIC\n")
	    gic <- gic.DMR(model.1percent)

	  } else  if (model_choice=="cv.DMRnet") {
	      cat("DMRnet with cv\n")
	      model.1percent <- tryCatch(cv_DMRnet(adult.train.1percent.x, adult.train.1percent.y, nlambda=100, family="binomial", nfolds=5, agressive = FALSE),
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
	  } else if (model_choice=="gic.GLAF") {
	    cat("GLAF method\n")
	    model.1percent <- tryCatch(glaf_4glm(adult.train.1percent.x, adult.train.1percent.y, nlambda=100),
	                               error=function(cond) {
	                                 message("Numerical instability in GLAF detected. Will skip this 1-percent set. Original error:")
	                                 message(cond)
	                                 return(list("red_light"))
	                               })

	    if (model.1percent[[1]] == "red_light") {
	      next
	    }

	    cat("GIC\n")
	    gic <- gic.DMR(model.1percent)   #we are using existing gic calculation which is compatible with GLAF models

	  } else  if (model_choice=="cv.GLAF") {
	    cat("GLAF with cv\n")
	    model.1percent <- tryCatch(cv_DMRnet(adult.train.1percent.x, adult.train.1percent.y, method="GLAF", nlambda=100, family="binomial", nfolds=5, agressive = FALSE),
	                               error=function(cond) {
	                                 message("Numerical instability in cv.GLAF detected. Will skip this 1-percent set. Original error:")
	                                 message(cond)
	                                 return(list("red_light"))
	                               })

	    if (model.1percent[[1]] == "red_light") {
	      next
	    }

	    #plot(model.1percent)
	    #gic <- gic.DMR(model.1percent, c = 2)
	    #plot(gic)
	  }  else if (model_choice=="scope") {
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




	  if (model_choice=="gic.DMRnet" | model_choice=="gic.GLAF") {
	    cat(model_choice, "pred\n")
	    prediction<- tryCatch(predict(model.1percent, newx=adult.test.1percent.x, df = gic$df.min, type="class"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (DMRnet) detected. Will skip this 1-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })

	    if (length(prediction)==2) {
	      next
	    }
	  } else  if (model_choice=="cv.DMRnet" | model_choice=="cv.GLAF") {
	    cat(model_choice, "pred\n")
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
	    prediction<- tryCatch(predict(model.1percent, newx=makeX(adult.test.1percent.x), s="lambda.min", type="class"),
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

	  if (model_choice == "gic.DMRnet" | model_choice == "gic.GLAF")
	    dfmin[run]<-gic$df.min
	  if (model_choice == "cv.DMRnet"  | model_choice == "cv.GLAF")
	    dfmin[run]<-model.1percent$df.min
	  if (model_choice == "cv.glmnet" )
	    dfmin[run]<-sum(coef(model.1percent, s="lambda.min")!=0)-1
	  if (model_choice == "scope")
	    dfmin[run]<-sum(abs(model.1percent$beta.best[[1]]) > 1e-10) +
	                sum(sapply(sapply(sapply(lapply(model.1percent$beta.best[[2]], as.factor), levels), unique), length)-1)
              	  #  length(unique(c(sapply(sapply(model.10percent$beta.best[[2]], as.factor), levels), sapply(sapply(model.10percent$beta.best[[1]], as.factor), levels),recursive=TRUE)))-1 + #-1 is for "0" level
              	  #             -sum(sapply(sapply(model.10percent$beta.best[[2]], as.factor), levels)!="0")   #and we subtract the number of factors = number of constraints from eq. (8) in Stokell et al.
              	  #the commented formula above had problems with levels close to 0 but nonzero, like these:

              	  #[[91]]
              	  #0                    1
              	  #6.28837260041593e-18 6.28837260041593e-18
              	  #Levels: 6.28837260041593e-18

              	  #[[92]]
              	  #0                    1
              	  #6.28837260041593e-18 6.28837260041593e-18
              	  #Levels: 6.28837260041593e-18
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
dev.off()

pdf("adult_errors.pdf",width=12,height=5)
boxplot(errors, ylim=c(0.16, 0.26))
dev.off()

pdf("adult_model_sizes.pdf",width=9,height=5)
boxplot(sizes)
dev.off()

pdf("adult_effective_lengths.pdf",width=12,height=5)
boxplot(effective_lengths)
dev.off()


