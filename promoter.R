

library(randomForest)
library(glmnet)
library(stats)  #glm
library(CatReg)
library(DMRnet)
library(digest)

source("glaf_4glm.R")
source("cv_DMRnet.R")

#library(devtools)
#load_all()

set.seed(strtoi(substr(digest("promoter", "md5", serialize = FALSE),1,7),16))



promoter.all<-read.csv("promoters.data", header=FALSE)
colnames <- c('class', "name", "nucleotides") #
colnames(promoter.all) <- colnames
promoter.all[,3] <- gsub("\\t", "", promoter.all[,3])  #removing tabs
promoter.all <- cbind(promoter.all, matrix(unlist(strsplit(promoter.all[,3], "")), ncol=57, byrow = TRUE))  #splitting nucleotides

for (i in c(1,4:60))
  promoter.all[,i]<-factor(promoter.all[,i])   #creating factors for response and for 57 nucleotide columns

levels(promoter.all[,1])[1]<-0
levels(promoter.all[,1])[2]<-1

cat("data loaded\n")

errors<-list()
effective_lengths<-list()
sizes<-list()
computation_times<-list()

gamma<-250

#1 PERCENT TRAIN / 99 PERCENT TEST SPLIT
runs<-1000
for (model_choice in c( "scope",  "cv.GLAF", "gic.GLAF", "cv.DMRnet", "gic.DMRnet", "cv.glmnet", "scope", "scope", "RF", "lr")) {
	gamma <- 350 - gamma    #it alternates between 250 and 100
	times<-dfmin<-misclassification_error<-lengths<-rep(0,runs)
	run<-1

	while (run<=runs) {
	  cat("generating train/test sets\n")
	  sample.70percent <- sample(1:nrow(promoter.all), 0.7*nrow(promoter.all))
	  promoter.train.70percent.x <- promoter.all[sample.70percent,4:60]
	  promoter.train.70percent.y <- promoter.all[sample.70percent,1]

	  promoter.test.70percent.x <- promoter.all[-sample.70percent,4:60]
	  promoter.test.70percent.y <- promoter.all[-sample.70percent,1]


	  #####RECOMPUTATION OF RELEVANT FACTORS in train set, to remove levels with no representative data (empty factors). Needed for random forest and glmnet
	  ###and for DMRnet - old package
	  ###but nor for DMRnet - new package
	  for (i in 1:57)
	    promoter.train.70percent.x[,i] <- factor(promoter.train.70percent.x[,i])

	  #remove data from test set with factors not present in train subsample as this causes predict() to fail
	  for (i in 1:57) {
	    train.levels <- levels(promoter.train.70percent.x[,i])
	    promoter.test.70percent.y<-promoter.test.70percent.y[which(promoter.test.70percent.x[,i] %in% train.levels)]
	    promoter.test.70percent.x<-promoter.test.70percent.x[which(promoter.test.70percent.x[,i] %in% train.levels),]
	  }
	  for (i in 1:57)
	    promoter.test.70percent.x[,i] <- factor(promoter.test.70percent.x[,i])   #recalculate factors now for new test


	  #removing columns with only one level:
	  singular_factors<-which(sapply(lapply(promoter.train.70percent.x, levels), length)==1)
	  if (length(singular_factors)>0) {
  	  promoter.test.70percent.x <- promoter.test.70percent.x[,-singular_factors]
  	  promoter.train.70percent.x <- promoter.train.70percent.x[,-singular_factors]
  	  cat("removed", length(singular_factors), "columns due to singular factors\n")
	  }

	  start.time <- Sys.time()
	  cat("Started: ", start.time,"\n")

	  if (model_choice=="gic.DMRnet") {
	    cat("DMRnet with GIC only\n")
	    model.70percent <- tryCatch(DMRnet(promoter.train.70percent.x, promoter.train.70percent.y, nlambda=100, family="binomial"),
	                               error=function(cond) {
	                                 message("Numerical instability in DMRnet detected. Will skip this 1-percent set. Original error:")
	                                 message(cond)
	                                 return(list("red_light"))
	                               })

	    if (model.70percent[[1]] == "red_light") {
	      next
	    }

	    cat("GIC\n")
	    gic <- gic.DMR(model.70percent, c = 2)

	  } else  if (model_choice=="cv.DMRnet") {
	      cat("DMRnet with cv\n")
	      model.70percent <- tryCatch(cv_DMRnet(promoter.train.70percent.x, promoter.train.70percent.y, nlambda=100, family="binomial", nfolds=5, agressive=FALSE),
	                                error=function(cond) {
	                                  message("Numerical instability in cv.DMRnet detected. Will skip this 1-percent set. Original error:")
	                                  message(cond)
	                                  return(list("red_light"))
	                                })

	    if (model.70percent[[1]] == "red_light") {
	      next
	    }

	    #plot(model.70percent)
	    #gic <- gic.DMR(model.70percent, c = 2)
	    #plot(gic)
	  } else if (model_choice=="gic.GLAF") {
	    cat("GLAF method\n")
	    model.70percent <- tryCatch(glaf_4glm(promoter.train.70percent.x, promoter.train.70percent.y, nlambda=100),
	                               error=function(cond) {
	                                 message("Numerical instability in GLAF detected. Will skip this 1-percent set. Original error:")
	                                 message(cond)
	                                 return(list("red_light"))
	                               })

	    if (model.70percent[[1]] == "red_light") {
	      next
	    }

	    cat("GIC\n")
	    gic <- gic.DMR(model.70percent, c = 2)   #we are using existing gic calculation which is compatible with GLAF models

	  } else  if (model_choice=="cv.GLAF") {
	    cat("GLAF with cv\n")
	    model.70percent <- tryCatch(cv_DMRnet(promoter.train.70percent.x, promoter.train.70percent.y, method="GLAF", nlambda=100, family="binomial", nfolds=5, agressive=FALSE),
	                               error=function(cond) {
	                                 message("Numerical instability in cv.DMRnet detected. Will skip this 1-percent set. Original error:")
	                                 message(cond)
	                                 return(list("red_light"))
	                               })

	    if (model.70percent[[1]] == "red_light") {
	      next
	    }

	    #plot(model.70percent)
	    #gic <- gic.DMR(model.70percent, c = 2)
	    #plot(gic)
	  }  else if (model_choice=="scope") {
	    cat("Scope, no cv, gamma=", gamma,"\n")
	    model.70percent <- tryCatch(scope.logistic(promoter.train.70percent.x, as.numeric(levels(promoter.train.70percent.y))[promoter.train.70percent.y], gamma=gamma),
	                               error=function(cond) {
	                                 message("Numerical instability in SCOPE detected. Will skip this 1-percent set. Original error:")
	                                 message(cond)
	                                 return(list("red_light"))
	                               })

	    if (model.70percent[[1]] == "red_light") {
	      next
	    }

	  } else if (model_choice=="RF") {
	    cat("random forest. no cv\n")
	    model.70percent <- randomForest(promoter.train.70percent.x, y=promoter.train.70percent.y)
	  } else if (model_choice=="lr") {
	    cat("Linear Regression no cv\n")
	    model.70percent <- glm(promoter.train.70percent.y~., data = promoter.train.70percent.x, family="binomial")
	  } else if (model_choice=="cv.glmnet") {
	    cat("glmnet with cv\n")
	    model.70percent<-cv.glmnet(makeX(promoter.train.70percent.x), promoter.train.70percent.y, family="binomial", nfolds=5)
	  } else
	    stop("Uknown method")




	  if (model_choice=="gic.DMRnet" | model_choice=="gic.GLAF") {
	    cat(model_choice, "pred\n")
	    prediction<- tryCatch(predict(model.70percent, newx=promoter.test.70percent.x, df = gic$df.min, type="class"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (DMRnet) detected. Will skip this 1-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })

	    if (length(prediction)==2) {
	      next
	    }
	  } else  if (model_choice=="cv.DMRnet" | model_choice =="cv.GLAF") {
	    cat(model_choice, "pred\n")
	    prediction<- tryCatch(predict(model.70percent, newx=promoter.test.70percent.x, type="class"),#df = gic$df.min, type="class"),
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
	    prediction<- ifelse(predict(model.70percent, promoter.test.70percent.x) >0.5,1,0)
	  } else if (model_choice=="RF") {
	    cat("Random Forest pred\n")
	    prediction<- tryCatch(predict(model.70percent, promoter.test.70percent.x, type="class"),
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
	    prediction<- ifelse(predict(model.70percent, promoter.test.70percent.x) >0,1,0)
	  } else if (model_choice=="cv.glmnet") {
	    cat("glmnet pred\n")
	    prediction<- tryCatch(predict(model.70percent, newx=makeX(promoter.test.70percent.x), s="lambda.min", type="class"),
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
	  misclassification_error[run]<-mean(prediction[!is.na(prediction)] != promoter.test.70percent.y[!is.na(prediction)])

	  if (model_choice == "gic.DMRnet" | model_choice == "gic.GLAF")
	    dfmin[run]<-gic$df.min
	  if (model_choice == "cv.DMRnet" | model_choice == "cv.GLAF")
	    dfmin[run]<-model.70percent$df.min
	  if (model_choice == "cv.glmnet" )
	    dfmin[run]<-sum(coef(model.70percent, s="lambda.min")!=0)-1
	  if (model_choice == "scope")
	    dfmin[run]<-sum(abs(model.70percent$beta.best[[1]]) > 1e-10) +
	                sum(sapply(sapply(sapply(lapply(model.70percent$beta.best[[2]], as.factor), levels), unique), length)-1)
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

write.csv(errors, "promoter_errors.csv")
write.csv(effective_lengths, "promoter_effective_lengths.csv")
write.csv(sizes, "promoter_model_sizes.csv")
write.csv(computation_times, "promoter_computation_times.csv")


pdf("promoter_computation_times.pdf",width=14,height=5)
boxplot(computation_times)
dev.off()

pdf("promoter_errors.pdf",width=14,height=5)
boxplot(errors)
dev.off()

pdf("promoter_model_sizes.pdf",width=12,height=5)
boxplot(sizes)
dev.off()

pdf("promoter_effective_lengths.pdf",width=14,height=5)
boxplot(effective_lengths)
dev.off()


