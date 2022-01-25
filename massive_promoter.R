


library(devtools)
library(vioplot)
load_all()


data("promoter")
promoter.errors<-read.csv("promoter_data/promoter_errors.csv")
promoter.df<-read.csv("promoter_data/promoter_model_sizes.csv")

cat("promoter data loaded\n")

runs<-200

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



    #removing columns with only one level, this would cause errors in DMRnet by design
    singular_factors<-which(sapply(lapply(promoter.train.70percent.x, levels), length)==1)
    if (length(singular_factors)>0) {
      promoter.test.70percent.x <- promoter.test.70percent.x[,-singular_factors]
      promoter.train.70percent.x <- promoter.train.70percent.x[,-singular_factors]
      cat("removed", length(singular_factors), "columns due to singular factors\n")
    }



    if (model_choice=="gic.DMRnet") {
      index=5
      if (run==1) {
        cat("DMRnet with GIC only\n")
      }
      model.70percent <- DMRnet(promoter.train.70percent.x, promoter.train.70percent.y, family="binomial")

      gic <- gic.DMR(model.70percent, c = 2)

    } else  if (model_choice=="cvg.DMRnet") {
      index=4
      if (run==1) {
        cat(model_choice, "with cvg\n")
      }
      model.70percent <- cv.DMRnet(promoter.train.70percent.x, promoter.train.70percent.y, family="binomial", indexation.mode = "GIC")

    } else if (model_choice=="gic.GLAMER") {
      index=3
      if (run==1) {
        cat("GLAMER method\n")
      }
      model.70percent <- DMRnet(promoter.train.70percent.x, promoter.train.70percent.y, algorithm="glamer", family="binomial", )

      gic <- gic.DMR(model.70percent, c = 2)   #we are using existing gic calculation which is compatible with GLAMER models

    }  else  if (model_choice=="cv+sd.GLAMER") {
      index=2
      if (run==1) {
        cat("GLAMER with cv+sd\n")
      }
      model.70percent<-cv.DMRnet(promoter.train.70percent.x, promoter.train.70percent.y, family="binomial", indexation.mode = "size", algorithm="glamer")

    } else
      stop("Uknown method")

    if (model_choice=="gic.DMRnet" | model_choice=="gic.GLAMER") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model.70percent, newx=promoter.test.70percent.x, df = gic$df.min, type="class")
      prediction1<- predict(gic, newx=promoter.test.70percent.x, type="class")
      cat("diff: ", sum(prediction1-prediction), "\n")
    } else  if (model_choice=="cvg.DMRnet") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model.70percent, newx=promoter.test.70percent.x, type="class")
    } else  if (model_choice=="cv+sd.GLAMER") {
      #cat(model_choice, "pred\n")
      prediction<- predict(model.70percent, newx=promoter.test.70percent.x, type="class", size="df.1se")
    } else
      stop("Uknown method")





    prediction[is.na(prediction)] <- 0
    misclassification_error<-mean(prediction[!is.na(prediction)] != promoter.test.70percent.y[!is.na(prediction)])

     if (model_choice == "gic.DMRnet" | model_choice == "gic.GLAMER")
      df<-gic$df.min
    if (model_choice == "cvg.DMRnet" )
      df<-model.70percent$df.min
    if (model_choice == "cv+sd.GLAMER" )
      df<-model.70percent$df.1se

    cat(run, "error = ", misclassification_error, "->", ecdf(promoter.errors[[index]])(misclassification_error), "\n")
    cat(run, "df = ", df, "->", ecdf(promoter.df[[index]])(df), "\n")

    mes[run]<-misclassification_error
    dfs[run]<-df

    run<-run+1
  }


  vioplot(list(mes, promoter.errors[[index]]), xlab = model_choice, ylab="error", main="promoter")
  vioplot(list(dfs, promoter.df[[index]]), xlab = model_choice, ylab="model size", main="promoter")
}

#write.csv(errors, "results/promoter_errors.csv")
#write.csv(sizes, "results/promoter_model_sizes.csv")

