
library(randomForest)
library(glmnet)
library(stats)  #model.matrix
library(CatReg)






insurance.all<-read.csv("train.csv", header=TRUE, comment.char="|", stringsAsFactors = TRUE)
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

for (model_choice in c( run_list )) {
	gamma <- 40 - gamma    #it alternates between 32 and 8
	times<-dfmin<-MSPE<-lengths<-rep(0,runs)
	run<-1

	while (run<=runs) {
	  cat("original response vector\n")

	  ################ original response 8 levels!
	  insurance.all.y<-insurance.all[,ncol(insurance.all)] +0.0
	  insurance.all.y.no_error<-insurance.all[,ncol(insurance.all)] +0.0

	  cat("generating train/test sets\n")

	  sample.10percent <- sample(1:nrow(insurance.all.x), 0.1*nrow(insurance.all.x))
    insurance.train.10percent.x <- insurance.all.x[sample.10percent,]
	  insurance.train.10percent.y <- insurance.all.y[sample.10percent]
	  insurance.test.10percent.x <- insurance.all.x[-sample.10percent,]
	  insurance.test.10percent.y <- insurance.all.y[-sample.10percent]
	  insurance.test.10percent.y.no_error <- insurance.all.y.no_error[-sample.10percent]


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
	      insurance.test.10percent.y.no_error<-insurance.test.10percent.y.no_error[which(insurance.test.10percent.x[,i] %in% train.levels)]
	      insurance.test.10percent.x<-insurance.test.10percent.x[which(insurance.test.10percent.x[,i] %in% train.levels),]
	    }

	  #recomputation of test factors - they should match train factors even if relevant levels not present
	  for (i in 1:ncol(insurance.train.10percent.x))
	    if (!(i %in% cont_columns)) {
	      train.levels <- levels(insurance.train.10percent.x[,i])
	      insurance.test.10percent.x[,i] <- factor(insurance.test.10percent.x[,i], levels=train.levels)   #recalculate factors now for new test
	    }

	  #HANDLING THE SINGULAR CASE

	  #removing columns with only one level:
	  singular_factors<-which(sapply(sapply(insurance.train.10percent.x, levels), length)==1)   #for continous columns length is 0
	  if (length(singular_factors)>0) {
  	  insurance.test.10percent.x <- insurance.test.10percent.x[,-singular_factors]
  	  insurance.train.10percent.x <- insurance.train.10percent.x[,-singular_factors]
  	  cat("removed", length(singular_factors), "columns due to singular factors\n")
	  }

	  insurance.make <- makeX(insurance.train.10percent.x)
	  prev_pos <- length(levels(insurance.train.10percent.x[,1]))
  	for (i in 2:ncol(insurance.train.10percent.x))   #we use the fact that we KNOW the first column is a "Product_Info_1" factor
	    if (!(i %in% cont_columns)) {  #removing columns from the last level, it is linearly dependent with the exception of the first column
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
	        if (i==1) {
	          reverse_lookup[pos:(pos+length(levels(insurance.train.10percent.x[,i]))-1)]<-i  #there are levels columns corresponding to the first original column
	          pos<-pos+length(levels(insurance.train.10percent.x[,i]))
	        } else {
	          reverse_lookup[pos:(pos+length(levels(insurance.train.10percent.x[,i]))-2)]<-i  #there are levels-1 columns corresponding to each original column other than the first
	          pos<-pos+length(levels(insurance.train.10percent.x[,i]))-1
	        }
	      } else {
	        reverse_lookup[pos]<-i
	        pos<-pos+1
	      }
	    #removal of columns for pivot positions larger than rank
	    remove_us<-reverse_lookup[QR$pivot[(QR$rank+1):length(QR$pivot)]]
	    insurance.train.10percent.x <- insurance.train.10percent.x[,-remove_us]
	    insurance.test.10percent.x <- insurance.test.10percent.x[,-remove_us]
	    cat("removed", length(unique(remove_us)), "columns (data linearly dependent)\n")
	  }

	  cat("consolidated factors and columns\n")

	  start.time <- Sys.time()
	  cat("Started: ", start.time,"\n")

	  if (model_choice=="cv.grLasso" | model_choice=="cv.MCP") {
	    cat(model_choice, "with CV\n")
	    X<-stats::model.matrix(~., insurance.train.10percent.x)
	    level_count <- sapply(lapply(insurance.train.10percent.x, levels), length)
	    level_count[level_count == 0] <- 2   #make it two for continous variables
	    groups<-rep(1:length(level_count), level_count-1)
	    if (model_choice == "cv.grLasso") {
	      penalty <-  "grLasso"
	    } else
	      penalty <- "grMCP"
	    model.10percent <- tryCatch(cv.grpreg(X[,-1], insurance.train.10percent.y, group=groups, penalty=penalty, nfolds=10),
	                                error=function(cond) {
	                                  message("Numerical instability in cv.grpreg detected. Will skip this 10-percent set. Original error:")
	                                  message(cond)
	                                  return(list("red_light"))
	                                })

	    if (model.10percent[[1]] == "red_light") {
	      next
	    }
	  } else if (model_choice=="cv.MCP-g") {
	    cat(model_choice, "with CV\n")
	    X<-stats::model.matrix(~., insurance.train.10percent.x)
	    level_count <- sapply(lapply(insurance.train.10percent.x, levels), length)
	    level_count[level_count == 0] <- 2   #make it two for continous variables
	    groups<-rep(1:length(level_count), level_count-1)
	    if (model_choice == "cv.grLasso") {
	      penalty <-  "grLasso"
	    } else
	      penalty <- "grMCP"
	    model.10percent <- tryCatch(cv.grpreg(X[,-1], insurance.train.10percent.y, group=groups, penalty=penalty, gamma=gamma, nfolds=10),
	                                error=function(cond) {
	                                  message("Numerical instability in cv.grpreg detected. Will skip this 10-percent set. Original error:")
	                                  message(cond)
	                                  return(list("red_light"))
	                                })

	    if (model.10percent[[1]] == "red_light") {
	      next
	    }
	  } else if (model_choice=="gic.DMRnet") {
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

	    cat("GIC\n")
	    gic <- gic.DMR(model.10percent)
	  } else  if (model_choice=="cvg.DMRnet"| model_choice == "cvg(e+m).DMRnet") {
	    cat(model_choice, "with cvg\n")
	    if (model_choice == "cvg(e+m).DMRnet") {
	      plateau_resistant <- TRUE
	    } else
	      plateau_resistant <- FALSE
	    model.10percent <- tryCatch(cvg_DMRnet(insurance.train.10percent.x, insurance.train.10percent.y, nlambda=100, family="gaussian", nfolds=10, agressive=TRUE, plateau_resistant_CV = plateau_resistant),
	                                error=function(cond) {
	                                  message("Numerical instability in cvg.DMRnet detected. Will skip this 10-percent set. Original error:")
	                                  message(cond)
	                                  return(list("red_light"))
	                                })

	    if (model.10percent[[1]] == "red_light") {
	      next
	    }
	  } else  if (model_choice=="gic.GLAMER") {
	    cat("GLAMER with GIC only\n")
	    model.10percent <- tryCatch(glamer_4lm(insurance.train.10percent.x, insurance.train.10percent.y, nlambda=100),
	                                error=function(cond) {
	                                  message("Numerical instability in gic.GLAMER detected. Will skip this 10-percent set. Original error:")
	                                  message(cond)
	                                  return(list("red_light"))
	                                })

	    if (model.10percent[[1]] == "red_light") {
	      next
	    }

	    cat("GIC\n")
	    gic <- gic.DMR(model.10percent)
	  } else  if (model_choice=="cv+sd.GLAMER") {
	    cat("GLAMER with cv+sd\n")
	    model.10percent <- tryCatch(cv_sd_glamer(insurance.train.10percent.x, insurance.train.10percent.y, nlambda=100, family="gaussian", nfolds=10, agressive=TRUE),
	                                error=function(cond) {
	                                  message("Numerical instability in cv+sd.GLAMER detected. Will skip this 10-percent set. Original error:")
	                                  message(cond)
	                                  return(list("red_light"))
	                                })

	    if (model.10percent[[1]] == "red_light") {
	      next
	    }
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

	  } else if (model_choice=="lr") {
	    cat("Linear Regression no cv\n")
	    model.10percent <- glm(insurance.train.10percent.y~., data = insurance.train.10percent.x, family="gaussian")
	  } else if (model_choice=="cv.glmnet") {
	    cat("glmnet with cv\n")
	    model.10percent<-cv.glmnet(makeX(insurance.train.10percent.x), insurance.train.10percent.y, family="gaussian", nfolds=10)
	  } else
	    stop("Uknown method")




	  if (model_choice=="cv.grLasso" | model_choice=="cv.MCP" | model_choice=="cv.MCP-g") {
	    cat(model_choice, "with CV prediction\n")
	    X_test<-stats::model.matrix(~., insurance.test.10percent.x)
	    prediction<- tryCatch(predict(model.10percent, X_test[,-1]),
	                          error=function(cond) {
	                            message("Numerical instability in predict (grpreg) detected. Will skip this 10-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })

	    if (length(prediction)==2) {
	      next
	    }
	  } else if (model_choice=="gic.DMRnet" | model_choice=="gic.GLAMER") {
	    cat(model_choice, "pred\n")
	    prediction<- tryCatch(predict(model.10percent, newx=insurance.test.10percent.x, df = gic$df.min, type="response"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (DMRnet) detected. Will skip this 10-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })
	    #for GLAMER we use a predict from DMR which is compatible with GLAMER models
	    if (length(prediction)==2) {
	      next
	    }
	  } else  if (model_choice=="cvg.DMRnet" | model_choice=="cv+sd.GLAMER" | model_choice == "cv.GLAMER" | model_choice == "cv.DMRnet" | model_choice == "cvg(e+m).DMRnet" | model_choice=="cp+sd.GLAMER") {
	    cat(model_choice, "pred\n")
	    prediction<- tryCatch(predict(model.10percent, newx=insurance.test.10percent.x, type="response"),#df = gic$df.min, type="class"),
	                          error=function(cond) {
	                            message("Numerical instability in predict (DMRnet) detected. Will skip this 10-percent set. Original error:")
	                            message(cond)
	                            return(c(1,1))
	                          })
	    #for GLAMER we use a predict from DMR which is compatible with GLAMER models
	    if (length(prediction)==2) {
	      next
	    }
	  } else if (model_choice=="scope") {
	    cat("scope pred\n")
	    prediction<- predict(model.10percent, insurance.test.10percent.x)
	  } else if (model_choice=="lr") {
	    cat("Linear Regression pred\n")
	    prediction<- predict(model.10percent, insurance.test.10percent.x)
	  } else if (model_choice=="cv.glmnet") {
	    cat("glmnet pred\n")
	    prediction<- tryCatch(predict(model.10percent, newx=makeX(insurance.test.10percent.x), s="lambda.min"),
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

	  MSPE[run]<-mean((prediction[!is.na(prediction)] - insurance.test.10percent.y.no_error[!is.na(prediction)])^2)

	  if (model_choice=="cv.grLasso" | model_choice=="cv.MCP" | model_choice=="cv.MCP-g")
	    dfmin[run]<-sum(coef(model.10percent)!=0)
	  if (model_choice == "gic.DMRnet" | model_choice=="gic.GLAMER")
	    dfmin[run]<-gic$df.min
	  if (model_choice == "cvg.DMRnet" | model_choice=="cv+sd.GLAMER" | model_choice == "cv.GLAMER" | model_choice == "cv.DMRnet" | model_choice == "cvg(e+m).DMRnet" | model_choice=="cp+sd.GLAMER")
	    dfmin[run]<-model.10percent$df.min
	  if (model_choice == "cv.glmnet" )
	    dfmin[run]<-sum(coef(model.10percent, s="lambda.min")!=0)-1
	  if (model_choice == "scope")
	    dfmin[run]<-sum(abs(model.10percent$beta.best[[1]]) > 1e-10) +
	                sum(sapply(sapply(sapply(lapply(model.10percent$beta.best[[2]], as.factor), levels), unique), length)-1)
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

write.csv(errors, paste("results/", part_filename_and_number, "_errors.csv", sep=""))
write.csv(sizes, paste("results/", part_filename_and_number, "_model_sizes.csv", sep=""))



