
library(devtools)
library(vioplot)
library(digest)
load_all()

set.seed(strtoi(substr(digest("antigua", "md5", serialize = FALSE),1,7),16))   #making all that I can to reproduce previous results of version 0.2.0+glamer from summer 2021 (part of preparation for AAAI'22)

run_list <- c("cv+sd.GLAMER","gic.GLAMER",  "cvg.DMRnet", "gic.DMRnet")
antigua.expected.errors<-read.csv("data_antigua/antigua_errors.csv")
antigua.expected.df<-read.csv("data_antigua/antigua_model_sizes.csv")


library(DAAG)
data(antigua)
antigua[antigua[,6] == -9999,6] = NA
antigua <- na.omit(antigua)
antigua.all.x <- antigua[, -c(1,7)]

cat("antigua data loaded\n")

runs<-200

postscript("result_antigua.pdf", horizontal=TRUE,onefile=FALSE)
par(mfrow=c(2,4))

for (model_choice in c( run_list )) {

  mes<-dfs<-rep(0,runs)
  run<-1

  while (run<=runs) {

    ################ original response
    antigua.all.y<-antigua[,ncol(antigua)] +0.0


   # cat("generating train/test sets\n")

    sample.70percent <- sample(1:nrow(antigua.all.x), 0.7*nrow(antigua.all.x))
    antigua.train.70percent.x <- antigua.all.x[sample.70percent,]
    antigua.train.70percent.y <- antigua.all.y[sample.70percent]
    antigua.test.70percent.x <- antigua.all.x[-sample.70percent,]
    antigua.test.70percent.y <- antigua.all.y[-sample.70percent]





    #HANDLING THE SINGULAR CASE

    #removing columns with only one value:
    constant_columns<-which( apply(antigua.train.70percent.x, 2, function(x) length(unique(x))) == 1)
    if (length(constant_columns)>0) {
      antigua.test.70percent.x <- antigua.test.70percent.x[,-constant_columns]
      antigua.train.70percent.x <- antigua.train.70percent.x[,-constant_columns]
      cat("removed", length(constant_columns), "columns due to constant values\n")
    }



    if (model_choice=="gic.DMRnet") {
      index=5
      if (run==1) {
        cat("DMRnet with GIC only\n")
      }
      model <- DMRnet(antigua.train.70percent.x, antigua.train.70percent.y, family="gaussian")

      gic <- gic.DMR(model, c = 2)

    } else  if (model_choice=="cvg.DMRnet") {
      index=4
      if (run==1) {
        cat(model_choice, "with cvg\n")
      }
      model <- cv.DMRnet(antigua.train.70percent.x, antigua.train.70percent.y, family="gaussian", indexation.mode = "GIC")

    } else if (model_choice=="gic.GLAMER") {
      index=3
      if (run==1) {
        cat("GLAMER method\n")
      }
      model <- DMRnet(antigua.train.70percent.x, antigua.train.70percent.y, algorithm="glamer", family="gaussian")

      gic <- gic.DMR(model, c = 2)   #we are using existing gic calculation which is compatible with GLAMER models

    }  else  if (model_choice=="cv+sd.GLAMER") {
      index=2
      if (run==1) {
        cat("GLAMER with cv+sd\n")
      }
      model<-cv.DMRnet(antigua.train.70percent.x, antigua.train.70percent.y, family="gaussian", indexation.mode = "dimension", algorithm="glamer")

    } else
      stop("Uknown method")

    if (model_choice=="gic.DMRnet" | model_choice=="gic.GLAMER") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model, newx=antigua.test.70percent.x, df = gic$df.min, type="class")
      prediction1<- predict(gic, newx=antigua.test.70percent.x, type="class")
      df<-gic$df.min
      cat("diff: ", sum(prediction1-prediction), "\n")
    } else  if (model_choice=="cvg.DMRnet") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model, newx=antigua.test.70percent.x, type="class")
      df<-model$df.min
    } else  if (model_choice=="cv+sd.GLAMER") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model, newx=antigua.test.70percent.x, type="class", md="df.1se")
      df<-model$df.1se
    } else
      stop("Uknown method")



   # prediction[is.na(prediction)] <- 0
    MSPE<-mean((prediction[!is.na(prediction)] - antigua.test.70percent.y[!is.na(prediction)])^2)

    cat(run, "error = ", MSPE, "->", ecdf(antigua.expected.errors[[index]])(MSPE), "\n")
    cat(run, "df = ", df, "->", ecdf(antigua.expected.df[[index]])(df), "\n")

    mes[run]<-MSPE
    dfs[run]<-df

    run<-run+1
  }


  vioplot(list(actual=mes, expected=antigua.expected.errors[[index]]), ylab="Error", main=model_choice)
  vioplot(list(actual=dfs, expected=antigua.expected.df[[index]]), ylab="Model Dimension", main=model_choice)
}

graphics.off()
