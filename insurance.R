
library(randomForest)
library(glmnet)
library(stats)  #glm
library(CatReg)
library(DMRnet)
library(digest)
library(scriptName)

#
#library(devtools)
#load_all()


filename_and_dirs <- current_filename()
part_filename_and_number <- substr(filename_and_dirs, nchar(filename_and_dirs)-11, nchar(filename_and_dirs))
print(part_filename_and_number)
set.seed(strtoi(substr(digest(part_filename_and_number, "md5", serialize = FALSE),1,7),16))
cat("seed set as md5 hash of the following string: ", part_filename_and_number,"\n\n")
cat("search for the result files postfixed with the same string: ", part_filename_and_number,"\n\n")

source("cv_DMRnet.R")

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
runs<-25
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
	  insurance.all.y.no_error <- (insurance.all.y - m)/std
	  insurance.all.y <- insurance.all.y.no_error + rnorm(length(insurance.all.y), 0, 1) #response was then scaled to have unit variance, after which standard normal noise was added.

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
	  for (i in 1:ncol(insurance.train.10percent.x))
	    if (!(i %in% cont_columns)) {
	      train.levels <- levels(insurance.train.10percent.x[,i])
	      insurance.test.10percent.x[,i] <- factor(insurance.test.10percent.x[,i], levels=train.levels)   #recalculate factors now for new test
	    }

	  #handling the singular case
	  insurance.make <- makeX(insurance.train.10percent.x)
	  prev_pos <- length(levels(insurance.train.10percent.x[,1]))
  	for (i in 2:ncol(insurance.train.10percent.x))
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

	  MSPE[run]<-mean((prediction[!is.na(prediction)] - insurance.test.10percent.y.no_error[!is.na(prediction)])^2)

	  if (model_choice == "gic.DMRnet")
	    dfmin[run]<-gic$df.min
	  if (model_choice == "cv.DMRnet" )
	    dfmin[run]<-model.10percent$df.min
	  if (model_choice == "scope")
	    dfmin[run]<-sum(abs(model.10percent$beta.best[[1]]) > 1e-10) +
	                sum(sapply(sapply(sapply(sapply(model.10percent$beta.best[[2]], as.factor), levels), unique), length)-1)
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


	computation_times[[model_name]]<-times
	effective_lengths[[model_name]]<-lengths
	if (length(dfmin[dfmin>0])>0)
		sizes[[model_name]]<-dfmin
	errors[[model_name]]<-MSPE

}

write.csv(errors, paste(part_filename_and_number, "_errors.csv", sep=""))
write.csv(effective_lengths, paste(part_filename_and_number, "_effective_lengths.csv", sep=""))
write.csv(sizes, paste(part_filename_and_number, "_model_sizes.csv", sep=""))
write.csv(computation_times, paste(part_filename_and_number, "_computation_times.csv", sep=""))


pdf(paste(part_filename_and_number, "_computation_times.pdf", sep=""),width=12,height=5)
boxplot(computation_times)
dev.off()

pdf(paste(part_filename_and_number, "_errors.pdf", sep=""),width=12,height=5)
boxplot(errors)
dev.off()

pdf(paste(part_filename_and_number, "_model_sizes.pdf", sep=""),width=9,height=5)
boxplot(sizes)
dev.off()

pdf(paste(part_filename_and_number, "_effective_lengths.pdf", sep=""),width=12,height=5)
boxplot(effective_lengths)
dev.off()


