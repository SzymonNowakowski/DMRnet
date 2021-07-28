
library(devtools)
library(randomForest)
library(glmnet)
library(stats)  #glm
library(CatReg)
load_all()
cat("loaded\n")

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
adult.train[,14] = factor(adult.train[,14])
adult.train[,7] = factor(adult.train[,7])
adult.train[,2] = factor(adult.train[,2])


adult.test<-subset(adult.test, adult.test[,14] != " ?")
adult.test<-subset(adult.test, adult.test[,7] != " ?")
adult.test<-subset(adult.test, adult.test[,2] != " ?")
adult.test[,14] = factor(adult.test[,14])
adult.test[,7] = factor(adult.test[,7])
adult.test[,2] = factor(adult.test[,2])

#consiliation of different level names in train and test sets (they end with '.' in test set)
levels(adult.test[,15])[1]<-0
levels(adult.test[,15])[2]<-1
levels(adult.train[,15])[1]<-0
levels(adult.train[,15])[2]<-1


adult.all<-rbind(adult.train, adult.test)
####HURRAY. In total 45222 observations (train+test) as in Stokell's paper

adult.train.x = adult.train[,c(1,2,4,6:10,13:14)]  #I exclude only education_num and fnlwgt and capital_gain & capital_loss
adult.train.y = adult.train[,15]


#ORIGINAL TEST/TEST SPLIT FAILS:
#model.original <- DMRnet(adult.train.x, adult.train.y, nlambda=100, family="binomial") ###comment this to go further


#1 PERCENT TRAIN / 99 PERCENT TEST SPLIT
runs<-200
model_choice <- "scope"
gamma <- 250
dfmin<-misclassification_error<-lengths<-rep(0,runs)

for (run in 1:runs) {
  sample.1percent <- sample(1:nrow(adult.all), 0.01*nrow(adult.all))
  adult.train.1percent.x <- adult.all[sample.1percent,c(1,2,4,6:10,13:14)]
  adult.train.1percent.y <- adult.all[sample.1percent,15]

  adult.test.1percent.x <- adult.all[-sample.1percent,c(1,2,4,6:10,13:14)]
  adult.test.1percent.y <- adult.all[-sample.1percent,15]

  ######WITH NO RECOMPUTAION OF FACTORS, grpreg FAILS with NA in empty factors
  #MOD <- DMRnet(adult.train.1percent.x, adult.train.1percent.y, nlambda=100, family="binomial") ###comment this to go further
  ####because this fails in grpreg
  #Error: Missing data (NA's) detected in X.  You must eliminate missing data (e.g., by removing cases, removing features, or imputation) before passing X to grpreg

  #####RECOMPUTATION OF RELEVANT FACTORS in train set, to remove data in test set. Needed for random forest and glmnet
  ###but nor for DMRnet (not any more)
  for (i in c(2:8,10))
    adult.train.1percent.x[,i] <- factor(adult.train.1percent.x[,i])

  if (model_choice=="DMRnet") {
      cat("DMRnet with cv\n")
      model.1percent <- tryCatch(cv.DMRnet(adult.train.1percent.x, adult.train.1percent.y, nlambda=100, family="binomial", nfolds=5),
                                error=function(cond) {
                                  message("Numerical instability in DMRnet detected. Will skip this 1-percent set. Original error:")
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
    model.1percent <- scope.logistic(adult.train.1percent.x, as.numeric(levels(adult.train.1percent.y))[adult.train.1percent.y], gamma=gamma)
  } else if (model_choice=="rf") {
    cat("random forest. no cv\n")
    model.1percent <- randomForest(adult.train.1percent.x, y=adult.train.1percent.y)
  } else if (model_choice=="lr") {
    cat("Linear Regression no cv\n")
    model.1percent <- glm(adult.train.1percent.y~., data = adult.train.1percent.x, family="binomial")
  } else if (model_choice=="glmnet") {
    cat("glmnet with cv\n")
    model.1percent<-cv.glmnet(makeX(adult.train.1percent.x), adult.train.1percent.y, family="binomial", nfolds=5)
  } else
    stop("Uknown method")



  #remove data from test with factors  not present in train subsample as this causes predict() to fail
  for (i in c(2:8,10)) {
    train.levels <- levels(adult.train.1percent.x[,i])
    adult.test.1percent.y<-adult.test.1percent.y[which(adult.test.1percent.x[,i] %in% train.levels)]
    adult.test.1percent.x<-adult.test.1percent.x[which(adult.test.1percent.x[,i] %in% train.levels),]
  }
  for (i in c(2:8,10))
    adult.test.1percent.x[,i] <- factor(adult.test.1percent.x[,i])   #recalculate factors now for new test

  if (model_choice=="DMRnet") {
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
  } else if (model_choice=="rf") {
    cat("Random Forest pred\n")
    prediction<- predict(model.1percent, adult.test.1percent.x, type="class")
  } else if (model_choice=="lr") {
    cat("Linear Regression pred\n")
    prediction<- ifelse(predict(model.1percent, adult.test.1percent.x) >0,1,0)
  } else if (model_choice=="glmnet") {
    cat("glmnet pred\n")
    prediction<-predict(model.1percent, newx=makeX(adult.test.1percent.x), type="class")
  } else
    stop("Uknown method")

  misclassification_error[run]<-sum(prediction[!is.na(prediction)] != adult.test.1percent.y[!is.na(prediction)]) / length(adult.test.1percent.y[!is.na(prediction)])
  lengths[run]<-length(prediction[!is.na(prediction)])
  if (model_choice == "DMRnet")
    dfmin[run]<-model.1percent$df.min
  if (model_choice == "scope")
    dfmin[run]<-length(unique(c(sapply(sapply(model.1percent$beta.best[[2]], as.factor), levels), recursive=TRUE)))-1  #-1 is for "0" level

  cat(run, "median = ", median(misclassification_error[misclassification_error>0]), "\n")
  cat(run, "df.min = ", mean(dfmin[misclassification_error>0]), "\n")
  cat(run, "lengths = ", mean(lengths[misclassification_error>0]), "\n")
}
boxplot(misclassification_error[misclassification_error!=0])
cat("overall median = ", median(misclassification_error[misclassification_error!=0]), "\n")

pdf("dfmin.pdf",width=7,height=5)
boxplot(dfmin[dfmin>0])
dev.off()

pdf("dfmin_error.pdf",width=7,height=5)
plot(dfmin[dfmin>0], misclassification_error[misclassification_error>0])
dev.off()

