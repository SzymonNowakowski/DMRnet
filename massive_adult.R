

library(devtools)
library(vioplot)
library(digest)
load_all()

set.seed(strtoi(substr(digest("adult", "md5", serialize = FALSE),1,7),16))


adult.train<-read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data", header=FALSE, comment.char="|", stringsAsFactors = TRUE)
adult.test<-read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.test", header=FALSE, comment.char="|", stringsAsFactors = TRUE)

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

adult.expected.errors<-read.csv("data_adult/adult_errors.csv")
adult.expected.df<-read.csv("data_adult/adult_model_sizes.csv")


cat("adult data loaded\n")

#1 PERCENT TRAIN / 99 PERCENT TEST SPLIT
runs<-200

postscript("result_adult.ps", horizontal=TRUE,onefile=FALSE)
par(mfrow=c(2,4))

for (model_choice in c( "cv+sd.PDMR", "gic.PDMR", "cvg.DMRnet", "gic.DMRnet")) {
  mes<-dfs<-rep(0,runs)
  run<-1

  while (run<=runs) {
    #cat("generating train/test sets\n")
    sample.1percent <- sample(1:nrow(adult.all), 0.01*nrow(adult.all))
    adult.train.1percent.x <- adult.all[sample.1percent,c(1,2,4,6:10,13:14)] #I exclude education_num and fnlwgt and capital_gain & capital_loss
    adult.train.1percent.y <- adult.all[sample.1percent,15]

    adult.test.1percent.x <- adult.all[-sample.1percent,c(1,2,4,6:10,13:14)]
    adult.test.1percent.y <- adult.all[-sample.1percent,15]




    #removing columns with only one value:
    constant_columns<-which( apply(adult.train.1percent.x, 2, function(x) length(unique(x))) == 1)

    if (length(constant_columns)>0) {
      adult.test.1percent.x <- adult.test.1percent.x[,-constant_columns]
      adult.train.1percent.x <- adult.train.1percent.x[,-constant_columns]
      cat("removed", length(constant_columns), "columns due to constant values\n")
    }


    if (model_choice=="gic.DMRnet") {
      index=5
      if (run==1) {
        cat("DMRnet with GIC only\n")
      }
      model <- DMRnet(adult.train.1percent.x, adult.train.1percent.y, family="binomial")

      gic <- gic.DMR(model, c = 2)

    } else  if (model_choice=="cvg.DMRnet") {
      index=4
      if (run==1) {
        cat(model_choice, "with cvg\n")
      }
      model <- cv.DMRnet(adult.train.1percent.x, adult.train.1percent.y, family="binomial", indexation.mode = "GIC")

    } else if (model_choice=="gic.PDMR") {
      index=3
      if (run==1) {
        cat("PDMR method\n")
      }
      model <- DMRnet(adult.train.1percent.x, adult.train.1percent.y, algorithm="PDMR", family="binomial")

      gic <- gic.DMR(model, c = 2)   #we are using existing gic calculation which is compatible with PDMR models

    }  else  if (model_choice=="cv+sd.PDMR") {
      index=2
      if (run==1) {
        cat("PDMR with cv+sd\n")
      }
      model<-cv.DMRnet(adult.train.1percent.x, adult.train.1percent.y, family="binomial", indexation.mode = "dimension", algorithm="PDMR")

    } else
      stop("Uknown method")



    if (model_choice=="gic.DMRnet" | model_choice=="gic.PDMR") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model, newx=adult.test.1percent.x, df = gic$df.min, type="class", unknown.factor.levels="NA")
      prediction1<- predict(gic, newx=adult.test.1percent.x, type="class", unknown.factor.levels="NA")
      df<-gic$df.min
      cat("diff: ", sum(prediction1[!is.na(prediction1)]-prediction[!is.na(prediction)]), "\n")
    } else  if (model_choice=="cvg.DMRnet") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model, newx=adult.test.1percent.x, type="class", unknown.factor.levels="NA")
      df<-model$df.min
    } else  if (model_choice=="cv+sd.PDMR") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model, newx=adult.test.1percent.x, type="class", md="df.1se", unknown.factor.levels="NA")  #the very first test case of cv+sd.PDMR fails because of unknown factor levels, if not set to "NA"
      df<-model$df.1se
    } else
      stop("Uknown method")


    #prediction[is.na(prediction)] <- 0
    misclassification_error<-mean(prediction[!is.na(prediction)] != adult.test.1percent.y[!is.na(prediction)])


    cat(run, "error = ", misclassification_error, "->", ecdf(adult.expected.errors[[index]])(misclassification_error), "\n")
    cat(run, "df = ", df, "->", ecdf(adult.expected.df[[index]])(df), "\n")

    mes[run]<-misclassification_error
    dfs[run]<-df

    run<-run+1
  }


  vioplot(list(actual=mes, expected=adult.expected.errors[[index]]), ylab="Error", main=model_choice)
  vioplot(list(actual=dfs, expected=adult.expected.df[[index]]), ylab="Model Dimension", main=model_choice)
}

graphics.off()
