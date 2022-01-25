
library(devtools)
library(vioplot)
load_all()

run_list = c("cv+sd.GLAMER","gic.GLAMER",  "cvg.DMRnet", "gic.DMRnet")
antigua.errors<-read.csv("antigua_data/antigua_errors.csv")
antigua.df<-read.csv("antigua_data/antigua_model_sizes.csv")


library(DAAG)
data(antigua)
antigua[antigua[,6] == -9999,6] = NA
antigua <- na.omit(antigua)
antigua.all.x <- antigua[, -c(1,7)]
cont_columns<-c(3,5)

cat("antigua data loaded\n")

runs<-200


for (model_choice in c( run_list )) {

  mes<-dfs<-rep(0,runs)
  run<-1

  while (run<=runs) {

    ################ original response 8 levels!
    antigua.all.y<-antigua[,ncol(antigua)] +0.0
    antigua.all.y.no_error<-antigua.all.y

   # cat("generating train/test sets\n")

    sample.70percent <- sample(1:nrow(antigua.all.x), 0.7*nrow(antigua.all.x))
    antigua.train.70percent.x <- antigua.all.x[sample.70percent,]
    antigua.train.70percent.y <- antigua.all.y[sample.70percent]
    antigua.test.70percent.x <- antigua.all.x[-sample.70percent,]
    antigua.test.70percent.y <- antigua.all.y[-sample.70percent]
    antigua.test.70percent.y.no_error <- antigua.all.y.no_error[-sample.70percent]




    #HANDLING THE SINGULAR CASE

    #removing columns with only one level:
    singular_factors<-which(sapply(sapply(antigua.train.70percent.x, levels), length)==1)   #for continous columns length is 0
    if (length(singular_factors)>0) {
      antigua.test.70percent.x <- antigua.test.70percent.x[,-singular_factors]
      antigua.train.70percent.x <- antigua.train.70percent.x[,-singular_factors]
      cat("removed", length(singular_factors), "columns due to singular factors\n")
    }



    if (model_choice=="gic.DMRnet") {
      index=5
      if (run==1) {
        cat("DMRnet with GIC only\n")
      }
      model.70percent <- DMRnet(antigua.train.70percent.x, antigua.train.70percent.y, family="gaussian")

      gic <- gic.DMR(model.70percent, c = 2)

    } else  if (model_choice=="cvg.DMRnet") {
      index=4
      if (run==1) {
        cat(model_choice, "with cvg\n")
      }
      model.70percent <- cv.DMRnet(antigua.train.70percent.x, antigua.train.70percent.y, family="gaussian", indexation.mode = "GIC")

    } else if (model_choice=="gic.GLAMER") {
      index=3
      if (run==1) {
        cat("GLAMER method\n")
      }
      model.70percent <- DMRnet(antigua.train.70percent.x, antigua.train.70percent.y, algorithm="glamer", family="gaussian", )

      gic <- gic.DMR(model.70percent, c = 2)   #we are using existing gic calculation which is compatible with GLAMER models

    }  else  if (model_choice=="cv+sd.GLAMER") {
      index=2
      if (run==1) {
        cat("GLAMER with cv+sd\n")
      }
      model.70percent<-cv.DMRnet(antigua.train.70percent.x, antigua.train.70percent.y, family="gaussian", indexation.mode = "size", algorithm="glamer")

    } else
      stop("Uknown method")

    if (model_choice=="gic.DMRnet" | model_choice=="gic.GLAMER") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model.70percent, newx=antigua.test.70percent.x, df = gic$df.min, type="class")
      prediction1<- predict(gic, newx=antigua.test.70percent.x, type="class")
      cat("diff: ", sum(prediction1-prediction), "\n")
    } else  if (model_choice=="cvg.DMRnet") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model.70percent, newx=antigua.test.70percent.x, type="class")
    } else  if (model_choice=="cv+sd.GLAMER") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model.70percent, newx=antigua.test.70percent.x, type="class", size="df.1se")
    } else
      stop("Uknown method")





    prediction[is.na(prediction)] <- 0
    MSPE<-mean((prediction[!is.na(prediction)] - antigua.test.70percent.y.no_error[!is.na(prediction)])^2)

    if (model_choice == "gic.DMRnet" | model_choice == "gic.GLAMER")
      df<-gic$df.min
    if (model_choice == "cvg.DMRnet" )
      df<-model.70percent$df.min
    if (model_choice == "cv+sd.GLAMER" )
      df<-model.70percent$df.1se

    cat(run, "error = ", MSPE, "->", ecdf(antigua.errors[[index]])(MSPE), "\n")
    cat(run, "df = ", df, "->", ecdf(antigua.df[[index]])(df), "\n")

    mes[run]<-MSPE
    dfs[run]<-df

    run<-run+1
  }


  vioplot(list(mes, antigua.errors[[index]]), xlab = model_choice, ylab="error", main="antigua")
  vioplot(list(dfs, antigua.df[[index]]), xlab = model_choice, ylab="model size", main="antigua")
}
