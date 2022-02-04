


library(devtools)
library(vioplot)
library(digest)
load_all()

set.seed(strtoi(substr(digest("promoter", "md5", serialize = FALSE),1,7),16))   #making all that I can to reproduce previous results of version 0.2.0+glamer from summer 2021 (part of preparation for AAAI'22)

data("promoter")
promoter.expected.errors<-read.csv("data_promoter/promoter_errors.csv")
promoter.expected.df<-read.csv("data_promoter/promoter_model_sizes.csv")

cat("promoter data loaded\n")

runs<-200

postscript("promoter_result.pdf", horizontal=TRUE,onefile=FALSE)
par(mfrow=c(2,4))

for (model_choice in c("cv+sd.GLAMER",  "gic.GLAMER","cvg.DMRnet", "gic.DMRnet")) {
  mes<-dfs<-rep(0,runs)
  run<-1

  while (run<=runs) {
 #   cat("generating train/test sets\n")
    sample.70percent <- sample(1:nrow(promoter), 0.7*nrow(promoter))
    promoter.train.70percent.x <- promoter[sample.70percent,2:58]
    promoter.train.70percent.y <- promoter[sample.70percent,1]

    promoter.test.70percent.x <- promoter[-sample.70percent,2:58]
    promoter.test.70percent.y <- promoter[-sample.70percent,1]



    #removing columns with only one value, this would cause errors in DMRnet by design
    singular_columns<-which( apply(promoter.train.70percent.x, 2, function(x) length(unique(x))) == 1)
    if (length(singular_columns)>0) {
      promoter.test.70percent.x <- promoter.test.70percent.x[,-singular_columns]
      promoter.train.70percent.x <- promoter.train.70percent.x[,-singular_columns]
      cat("removed", length(singular_columns), "columns due to singular values\n")
    }



    if (model_choice=="gic.DMRnet") {
      index=5
      if (run==1) {
        cat("DMRnet with GIC only\n")
      }
      model <- DMRnet(promoter.train.70percent.x, promoter.train.70percent.y, family="binomial")

      gic <- gic.DMR(model, c = 2)

    } else  if (model_choice=="cvg.DMRnet") {
      index=4
      if (run==1) {
        cat(model_choice, "with cvg\n")
      }
      model <- cv.DMRnet(promoter.train.70percent.x, promoter.train.70percent.y, family="binomial", indexation.mode = "GIC")

    } else if (model_choice=="gic.GLAMER") {
      index=3
      if (run==1) {
        cat("GLAMER method\n")
      }
      model <- DMRnet(promoter.train.70percent.x, promoter.train.70percent.y, algorithm="glamer", family="binomial")

      gic <- gic.DMR(model, c = 2)   #we are using existing gic calculation which is compatible with GLAMER models

    }  else  if (model_choice=="cv+sd.GLAMER") {
      index=2
      if (run==1) {
        cat("GLAMER with cv+sd\n")
      }
      model<-cv.DMRnet(promoter.train.70percent.x, promoter.train.70percent.y, family="binomial", indexation.mode = "dimension", algorithm="glamer")

    } else
      stop("Uknown method")

    if (model_choice=="gic.DMRnet" | model_choice=="gic.GLAMER") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model, newx=promoter.test.70percent.x, df = gic$df.min, type="class")
      prediction1<- predict(gic, newx=promoter.test.70percent.x, type="class")
      df<-gic$df.min
      cat("diff: ", sum(prediction1-prediction), "\n")
    } else  if (model_choice=="cvg.DMRnet") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model, newx=promoter.test.70percent.x, type="class")
      df<-model$df.min
    } else  if (model_choice=="cv+sd.GLAMER") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model, newx=promoter.test.70percent.x, type="class", md="df.1se")
      df<-model$df.1se
    } else
      stop("Uknown method")



    #prediction[is.na(prediction)] <- 0
    misclassification_error<-mean(prediction[!is.na(prediction)] != promoter.test.70percent.y[!is.na(prediction)])

    cat(run, "error = ", misclassification_error, "->", ecdf(promoter.expected.errors[[index]])(misclassification_error), "\n")
    cat(run, "df = ", df, "->", ecdf(promoter.expected.df[[index]])(df), "\n")

    mes[run]<-misclassification_error
    dfs[run]<-df

    run<-run+1
  }


  vioplot(list(actual=mes, expected=promoter.expected.errors[[index]]), xlab = model_choice, ylab="error", main="promoter")
  vioplot(list(actual=dfs, expected=promoter.expected.df[[index]]), xlab = model_choice, ylab="model size", main="promoter")
}

graphics.off()

