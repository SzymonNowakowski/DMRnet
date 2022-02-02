
library(devtools)
library(vioplot)
library(digest)
load_all()

return_seed<-function(text, letter) { strtoi(substr(digest(paste("insurance_", text, "_orig", letter, ".R", sep=""), "md5", serialize = FALSE),1,7),16)}


run_list_glamer <- c("cv+sd.GLAMER","gic.GLAMER")
run_list_dmrnet <- c("cvg.DMRnet", "gic.DMRnet")
insurance.expected.errors<-read.csv("data_insurance/insurance_orig_total_errors.csv")
insurance.expected.df<-read.csv("data_insurance/insurance_orig_total_model_sizes.csv")



insurance.all<-read.csv("data_insurance/train.csv", header=TRUE, comment.char="|", stringsAsFactors = TRUE)
insurance.all<-insurance.all[,apply(apply(insurance.all,2,is.na), 2, sum)==0]  #removing columns with NA
insurance.all.x<-insurance.all[,2:(ncol(insurance.all)-1)]  #without ID and the response columns
cont_columns = c(4, 8, 9, 10, 11)

for (i in 1:ncol(insurance.all.x))
  if (!(i %in% cont_columns)) {
    insurance.all.x[,i] <- factor(insurance.all.x[,i])  #int->factors
  }

cat("insurance data loaded\n")

results_mes<-results_dfs<-matrix(nrow=200, ncol=4)

for (mega_run in 1:20) {
  for (part in c("glamer", "DMRnet")) {
    if (part == "glamer") {
      run_list<-run_list_glamer
    } else{
      run_lis<-run_list_dmrnet
    }

    set.seed(return_seed(part, letters[mega_run]))

    for (model_choice in c( run_list )) {


      for (run in 1:10) {

        ################ original response 8 levels!
        insurance.all.y<-insurance.all[,ncol(insurance.all)] +0.0

        cat("generating train/test sets\n")

        sample.10percent <- sample(1:nrow(insurance.all.x), 0.1*nrow(insurance.all.x))
        insurance.train.10percent.x <- insurance.all.x[sample.10percent,]
        insurance.train.10percent.y <- insurance.all.y[sample.10percent]
        insurance.test.10percent.x <- insurance.all.x[-sample.10percent,]
        insurance.test.10percent.y <- insurance.all.y[-sample.10percent]





        #HANDLING THE SINGULAR CASE

        #removing columns with only one value:
        singular_factors<-which( apply(insurance.train.10percent.x, 2, function(x) length(unique(x))) == 1)
        if (length(singular_factors)>0) {
          insurance.test.10percent.x <- insurance.test.10percent.x[,-singular_factors]
          insurance.train.10percent.x <- insurance.train.10percent.x[,-singular_factors]
          cat("removed", length(singular_factors), "columns due to singular values\n")
        }



        if (model_choice=="gic.DMRnet") {
          index=5
          if (run==1) {
            cat("DMRnet with GIC only\n")
          }
          model <- DMRnet(insurance.train.10percent.x, insurance.train.10percent.y, family="gaussian")

          gic <- gic.DMR(model, c = 2)

        } else  if (model_choice=="cvg.DMRnet") {
          index=4
          if (run==1) {
            cat(model_choice, "with cvg\n")
          }
          model <- cv.DMRnet(insurance.train.10percent.x, insurance.train.10percent.y, family="gaussian", indexation.mode = "GIC")

        } else if (model_choice=="gic.GLAMER") {
          index=3
          if (run==1) {
            cat("GLAMER method\n")
          }
          model <- DMRnet(insurance.train.10percent.x, insurance.train.10percent.y, algorithm="glamer", family="gaussian", )

          gic <- gic.DMR(model, c = 2)   #we are using existing gic calculation which is compatible with GLAMER models

        }  else  if (model_choice=="cv+sd.GLAMER") {
          index=2
          if (run==1) {
            cat("GLAMER with cv+sd\n")
          }
          model<-cv.DMRnet(insurance.train.10percent.x, insurance.train.10percent.y, family="gaussian", indexation.mode = "dimension", algorithm="glamer")

        } else
          stop("Uknown method")

        if (model_choice=="gic.DMRnet" | model_choice=="gic.GLAMER") {
          #cat(model_choice, "pred\n")
          prediction<- predict(model, newx=insurance.test.10percent.x, df = gic$df.min, type="class")
          prediction1<- predict(gic, newx=insurance.test.10percent.x, type="class")
          df<-gic$df.min
          cat("diff: ", sum(prediction1-prediction), "\n")
        } else  if (model_choice=="cvg.DMRnet") {
          #cat(model_choice, "pred\n")
          prediction<- predict(model, newx=insurance.test.10percent.x, type="class")
          df<-model$df.min
        } else  if (model_choice=="cv+sd.GLAMER") {
          #cat(model_choice, "pred\n")
          prediction<- predict(model, newx=insurance.test.10percent.x, type="class", md="df.1se")
          df<-model$df.1se
        } else
          stop("Uknown method")





        prediction[is.na(prediction)] <- 0
        MSPE<-mean((prediction[!is.na(prediction)] - insurance.test.10percent.y[!is.na(prediction)])^2)


        total_run<-run+10*(mega_run-1)
        cat(total_run, "error = ", MSPE, "->", ecdf(insurance.expected.errors[[index]])(MSPE), "\n")
        cat(total_run, "df = ", df, "->", ecdf(insurance.expected.df[[index]])(df), "\n")

        results_mes[total_run, index-1]<-MSPE
        results_dfs[total_run, index-1]<-df

      }#run
    }#model_choice
  }#part
}#mega_run


for (model_index in 1:4) {
  model_choice <- c("cv+sd.GLAMER","gic.GLAMER",  "cvg.DMRnet", "gic.DMRnet") [model_index]
  index <- model_index + 1 # a number of the corresponding expected result in a result filename
  vioplot(list(actual=results_mes[,model_index], expected=insurance.expected.errors[[index]]), xlab = model_choice, ylab="error", main="insurance")
  vioplot(list(actual=results_dfs[,model_index], expected=insurance.expected.df[[index]]), xlab = model_choice, ylab="model size", main="insurance")
}


