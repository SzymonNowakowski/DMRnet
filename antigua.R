
library(randomForest)
library(glmnet)
library(stats)  # model.matrix
library(CatReg)
library(DMRnet)
library(grpreg)
library(digest)


set.seed(strtoi(substr(digest("antigua", "md5", serialize = FALSE),1,7),16))

runs<-100
run_list = c(   "cv+sd.GLAMER","gic.GLAMER",  "cvg.DMRnet", "gic.DMRnet", "scope", "scope",  "cv.glmnet", "RF", "cv.MCP", "cv.MCP-g", "cv.MCP-g", "cv.grLasso")


source("cvg_DMRnet.R")
source("cv_sd_glamer.R")
source("glamer_4lm.R")

library(DAAG)
data(antigua)
antigua[antigua[,6] == -9999,6] = NA
antigua <- na.omit(antigua)
antigua.all.x <- antigua[, -c(1,7)]
cont_columns<-c(3,5)

cat("data loaded\n")

errors<-list()
effective_lengths<-list()
sizes<-list()
computation_times<-list()

gamma<-8


for (model_choice in c( run_list )) {
	gamma <- 40 - gamma    #it alternates between 32 and 8
	times<-dfmin<-MSPE<-lengths<-rep(0,runs)
	run<-1

	while (run<=runs) {
	  cat("original response vector\n")

	  ################ original response 8 levels!
	  antigua.all.y<-antigua[,ncol(antigua)] +0.0
	  antigua.all.y.no_error<-antigua.all.y

	  cat("generating train/test sets\n")

	  sample.70percent <- sample(1:nrow(antigua.all.x), 0.7*nrow(antigua.all.x))
    antigua.train.70percent.x <- antigua.all.x[sample.70percent,]
	  antigua.train.70percent.y <- antigua.all.y[sample.70percent]
	  antigua.test.70percent.x <- antigua.all.x[-sample.70percent,]
	  antigua.test.70percent.y <- antigua.all.y[-sample.70percent]
	  antigua.test.70percent.y.no_error <- antigua.all.y.no_error[-sample.70percent]


	  #####RECOMPUTATION OF RELEVANT FACTORS in train set, to remove levels with no representative data (empty factors). Needed for random forest and glmnet
	  ###and for DMRnet - old package
	  ###but nor for DMRnet - new package


	  for (i in 1:ncol(antigua.train.70percent.x))
	    if (!(i %in% cont_columns))
	      antigua.train.70percent.x[,i] <- factor(antigua.train.70percent.x[,i])


	  #remove data from test set with factors not present in train subsample as this causes predict() to fail
	  for (i in 1:ncol(antigua.train.70percent.x))
	    if (!(i %in% cont_columns)) {
	      train.levels <- levels(antigua.train.70percent.x[,i])
	      antigua.test.70percent.y<-antigua.test.70percent.y[which(antigua.test.70percent.x[,i] %in% train.levels)]
	      antigua.test.70percent.y.no_error<-antigua.test.70percent.y.no_error[which(antigua.test.70percent.x[,i] %in% train.levels)]
	      antigua.test.70percent.x<-antigua.test.70percent.x[which(antigua.test.70percent.x[,i] %in% train.levels),]
	    }

	  #recomputation of test factors - they should match train factors even if relevant levels not present
	  for (i in 1:ncol(antigua.train.70percent.x))
	    if (!(i %in% cont_columns)) {
	      train.levels <- levels(antigua.train.70percent.x[,i])
	      antigua.test.70percent.x[,i] <- factor(antigua.test.70percent.x[,i], levels=train.levels)   #recalculate factors now for new test
	    }

	  #HANDLING THE SINGULAR CASE

	  #removing columns with only one level:
	  singular_factors<-which(sapply(sapply(antigua.train.70percent.x, levels), length)==1)   #for continous columns length is 0
	  if (length(singular_factors)>0) {
  	  antigua.test.70percent.x <- antigua.test.70percent.x[,-singular_factors]
  	  antigua.train.70percent.x <- antigua.train.70percent.x[,-singular_factors]
  	  cat("removed", length(singular_factors), "columns due to singular factors\n")
	  }

	  cat("consolidated factors\n")

	  start.time <- Sys.time()
	  cat("Started: ", start.time,"\n")

	  if (model_choice=="cv.grLasso" | model_choice=="cv.MCP") {
	    cat(model_choice, "with CV\n")
	    X<-stats::model.matrix(~., antigua.train.70percent.x)
	    level_count <- sapply(lapply(antigua.train.70percent.x, levels), length)
	    level_count[level_count == 0] <- 2   #make it two for continous variables
	    groups<-rep(1:length(level_count), level_count-1)
	    if (model_choice == "cv.grLasso") {
	      penalty <-  "grLasso"
	    } else
	      penalty <- "grMCP"
	    model.70percent <- cv.grpreg(X[,-1], antigua.train.70percent.y, group=groups, penalty=penalty, nfolds=10)

	  } else if (model_choice=="cv.MCP-g") {
	    cat(model_choice, "with CV\n")
	    X<-stats::model.matrix(~., antigua.train.70percent.x)
	    level_count <- sapply(lapply(antigua.train.70percent.x, levels), length)
	    level_count[level_count == 0] <- 2   #make it two for continous variables
	    groups<-rep(1:length(level_count), level_count-1)
	    if (model_choice == "cv.grLasso") {
	      penalty <-  "grLasso"
	    } else
	      penalty <- "grMCP"
	    model.70percent <- cv.grpreg(X[,-1], antigua.train.70percent.y, group=groups, penalty=penalty, gamma=gamma, nfolds=10)

	  } else if (model_choice=="gic.DMRnet") {
	    cat("DMRnet with GIC only\n")
	    model.70percent <- tryCatch(DMRnet(antigua.train.70percent.x, antigua.train.70percent.y, nlambda=100, family="gaussian"),
	                               error=function(cond) {
	                                 message("Numerical instability in DMRnet detected. Will skip this 70-percent set. Original error:")
	                                 message(cond)
	                                 return(list("red_light"))
	                               })

	    if (model.70percent[[1]] == "red_light") {
	      next
	    }

	    cat("GIC\n")
	    gic <- gic.DMR(model.70percent)
	  } else  if (model_choice=="cvg.DMRnet" | model_choice == "cvg(e+m).DMRnet") {
	    cat(model_choice, "with cvg\n")
	    if (model_choice == "cvg(e+m).DMRnet") {
	      plateau_resistant <- TRUE
	    } else
	      plateau_resistant <- FALSE
	    model.70percent <- tryCatch(cvg_DMRnet(antigua.train.70percent.x, antigua.train.70percent.y, nlambda=100, family="gaussian", nfolds=10, plateau_resistant_CV =  plateau_resistant),
	                                error=function(cond) {
	                                  message("Numerical instability in cvg.DMRnet detected. Will skip this 70-percent set. Original error:")
	                                  message(cond)
	                                  return(list("red_light"))
	                                })

	    if (model.70percent[[1]] == "red_light") {
	      next
	    }
	  } else  if (model_choice=="gic.GLAMER") {
	    cat("GLAMER with GIC only\n")
	    model.70percent <- tryCatch(glamer_4lm(antigua.train.70percent.x, antigua.train.70percent.y, nlambda=100),
	                                error=function(cond) {
	                                  message("Numerical instability in gic.GLAMER detected. Will skip this 70-percent set. Original error:")
	                                  message(cond)
	                                  return(list("red_light"))
	                                })

	    if (model.70percent[[1]] == "red_light") {
	      next
	    }

	    cat("GIC\n")
	    gic <- gic.DMR(model.70percent)
	  } else  if (model_choice=="cv+sd.GLAMER") {
	    cat("GLAMER with cv+sd\n")

	    model.70percent <- tryCatch(cv_sd_glamer(antigua.train.70percent.x, antigua.train.70percent.y, nlambda=100, family="gaussian", nfolds=10),

	                                error=function(cond) {
	                                  message("Numerical instability in cv+sd.GLAMER detected. Will skip this 70-percent set. Original error:")
	                                  message(cond)
	                                  return(list("red_light"))
	                                })

	    if (model.70percent[[1]] == "red_light") {
	      next
	    }
	  } else if (model_choice=="scope") {
	    cat("Scope, no cv, gamma=", gamma,"\n")
	    model.70percent <- tryCatch(scope(antigua.train.70percent.x, antigua.train.70percent.y, gamma=gamma),
	                               error=function(cond) {
	                                 message("Numerical instability in SCOPE detected. Will skip this 70-percent set. Original error:")
	                                 message(cond)
	                                 return(list("red_light"))
	                               })

	    if (model.70percent[[1]] == "red_light") {
	      next
	    }

	  } else if (model_choice=="RF") {
	    cat("random forest. no cv\n")
	    model.70percent <- randomForest(antigua.train.70percent.x, y=antigua.train.70percent.y)
	  } else if (model_choice=="lr") {
	    cat("Linear Regression no cv\n")
	    model.70percent <- glm(antigua.train.70percent.y~., data = antigua.train.70percent.x, family="gaussian")
	  } else if (model_choice=="cv.glmnet") {
	    cat("glmnet with cv\n")
	    model.70percent<-cv.glmnet(makeX(antigua.train.70percent.x), antigua.train.70percent.y, family="gaussian", nfolds=10)
	  } else
	    stop("Uknown method")



	  if (model_choice=="cv.grLasso" | model_choice=="cv.MCP" | model_choice=="cv.MCP-g") {
	    cat(model_choice, "with CV prediction\n")
	    X_test<-stats::model.matrix(~., antigua.test.70percent.x)
	    prediction<-predict(model.70percent, X_test[,-1])
	  } else if (model_choice=="gic.DMRnet" | model_choice=="gic.GLAMER") {
	    cat(model_choice, "pred\n")
	    prediction<- tryCatch(predict(model.70percent, newx=antigua.test.70percent.x, df = gic$df.min, type="response"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (DMRnet) detected. Will skip this 70-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })
	    #for GLAMER we use a predict from DMR which is compatible with GLAMER models
	    if (length(prediction)==2) {
	      next
	    }
	  } else  if (model_choice=="cvg.DMRnet" | model_choice=="cv+sd.GLAMER" | model_choice == "cv.GLAMER" | model_choice == "cv.DMRnet" | model_choice == "cvg(e+m).DMRnet" | model_choice == "cp+sd.GLAMER") {
	    cat(model_choice, "pred\n")
	    prediction<- tryCatch(predict(model.70percent, newx=antigua.test.70percent.x, type="response"),#df = gic$df.min, type="class"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (DMRnet) detected. Will skip this 70-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })
	    #for GLAMER we use a predict from DMR which is compatible with GLAMER models
	    if (length(prediction)==2) {
	      next
	    }
	  } else if (model_choice=="scope") {
	    cat("scope pred\n")
	    prediction<- predict(model.70percent, antigua.test.70percent.x)
	  } else if (model_choice=="RF") {
	    cat("Random Forest pred\n")
	    prediction<- tryCatch(predict(model.70percent, antigua.test.70percent.x, type="class"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (RF) detected. Will skip this 70-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })

	    if (length(prediction)==2) {
	      next
	    }
	  } else if (model_choice=="lr") {
	    cat("Linear Regression pred\n")
	    prediction<- predict(model.70percent, antigua.test.70percent.x)
	  } else if (model_choice=="cv.glmnet") {
	    cat("glmnet pred\n")
	    prediction<- tryCatch(predict(model.70percent, newx=makeX(antigua.test.70percent.x), s="lambda.min"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (cv.glmnet) detected. Will skip this 70-percent set. Original error:")
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

	  MSPE[run]<-mean((prediction[!is.na(prediction)] - antigua.test.70percent.y.no_error[!is.na(prediction)])^2)

	  if (model_choice=="cv.grLasso" | model_choice=="cv.MCP" | model_choice=="cv.MCP-g")
	    dfmin[run]<-sum(coef(model.70percent)!=0)
	  if (model_choice == "gic.DMRnet" | model_choice=="gic.GLAMER")
	    dfmin[run]<-gic$df.min
	  if (model_choice == "cvg.DMRnet" | model_choice=="cv+sd.GLAMER" | model_choice == "cv.GLAMER" | model_choice == "cv.DMRnet" | model_choice == "cvg(e+m).DMRnet" | model_choice == "cp+sd.GLAMER")
	    dfmin[run]<-model.70percent$df.min
	  if (model_choice == "cv.glmnet" )
	    dfmin[run]<-sum(coef(model.70percent, s="lambda.min")!=0)-1
	  if (model_choice == "scope")
	    dfmin[run]<-sum(abs(model.70percent$beta.best[[1]]) > 1e-10) +
	                sum(sapply(sapply(sapply(lapply(model.70percent$beta.best[[2]], as.factor), levels), unique), length)-1)
	  #  length(unique(c(sapply(sapply(model.70percent$beta.best[[2]], as.factor), levels), sapply(sapply(model.70percent$beta.best[[1]], as.factor), levels),recursive=TRUE)))-1 + #-1 is for "0" level
	   #             -sum(sapply(sapply(model.70percent$beta.best[[2]], as.factor), levels)!="0")   #and we subtract the number of factors = number of constraints from eq. (8) in Stokell et al.
                       #the commented formula above had problems with levels close to 0 but nonzero, like these:

                  	  #[[91]]
                  	  #0                    1
                  	  #6.28837260041593e-18 6.28837260041593e-18
                  	  #Levels: 6.28837260041593e-18

                  	  #[[92]]
                  	  #0                    1
                  	  #6.28837260041593e-18 6.28837260041593e-18
                  	  #Levels: 6.28837260041593e-18

	  cat(run, "median = ", median(MSPE[MSPE>0]), "\n")
	  cat(run, "df.min = ", mean(dfmin[MSPE>0]), "\n")
	  cat(run, "lengths = ", mean(lengths[MSPE>0]), "\n")

	  run<-run+1
	}

	cat("overall median = ", median(MSPE[MSPE!=0]), "\n")


	model_name<-model_choice
	if (model_choice == "scope")
	  model_name<-paste(model_name, gamma, sep="-")

	if (model_choice == "cv.MCP-g")
	  model_name<-paste(model_name, gamma, sep="-")


	computation_times[[model_name]]<-times
	effective_lengths[[model_name]]<-lengths
	if (length(dfmin[dfmin>0])>0)
		sizes[[model_name]]<-dfmin
	errors[[model_name]]<-MSPE

}

write.csv(errors, "results/antigua_errors.csv")
write.csv(sizes, "results/antigua_model_sizes.csv")





