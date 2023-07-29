
library(devtools)
library(vioplot)
library(digest)
load_all()

return_seed<-function(text, letter) { strtoi(substr(digest(paste("insurance_", text, "_orig", letter, ".R", sep=""), "md5", serialize = FALSE),1,7),16)}


run_list_PDMR <- c("cv+sd.PDMR","gic.PDMR")
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

for (mega_run in 1:4) {
  for (part in c("PDMR", "DMRnet")) {
    if (part == "PDMR") {
      run_list<-run_list_PDMR
    } else{
      run_list<-run_list_dmrnet
    }

    set.seed(return_seed(part, letters[mega_run]))

    for (model_choice in c( run_list )) {


      for (run in 1:10) {

        ################ original response 8 levels!
        insurance.all.y<-insurance.all[,ncol(insurance.all)] +0.0

       #cat("generating train/test sets\n")

        sample.10percent <- sample(1:nrow(insurance.all.x), 0.1*nrow(insurance.all.x))
        insurance.train.10percent.x <- insurance.all.x[sample.10percent,]
        insurance.train.10percent.y <- insurance.all.y[sample.10percent]
        insurance.test.10percent.x <- insurance.all.x[-sample.10percent,]
        insurance.test.10percent.y <- insurance.all.y[-sample.10percent]





        #HANDLING THE SINGULAR CASE

        #removing columns with only one value:
        constant_columns<-which( apply(insurance.train.10percent.x, 2, function(x) length(unique(x))) == 1)
        if (length(constant_columns)>0) {
          insurance.test.10percent.x <- insurance.test.10percent.x[,-constant_columns]
          insurance.train.10percent.x <- insurance.train.10percent.x[,-constant_columns]
          cat("removed", length(constant_columns), "columns due to constant values\n")
        }



        if (model_choice=="gic.DMRnet") {
          index=5
          if (run==1) {
            cat("DMRnet with GIC only\n")
          }
          model <- tryCatch(DMRnet(insurance.train.10percent.x, insurance.train.10percent.y, family="gaussian"),
                            error=function(cond) {
                              message("Numerical instability in regular model creation detected. Will skip this set. Original error:")
                              message(cond); cat("\n")
                              return(c(1,1))
                            })
          if (length(model)==2) {
            next
          }

          gic <- gic.DMR(model, c = 2)

        } else  if (model_choice=="cvg.DMRnet") {
          index=4
          if (run==1) {
            cat(model_choice, "with cvg\n")
          }
          model <- tryCatch(cv.DMRnet(insurance.train.10percent.x, insurance.train.10percent.y, family="gaussian", indexation.mode = "GIC"),
                            error=function(cond) {
                              message("Numerical instability in model creation in CV (cv.DMRnet) detected. Will skip this set. Original error:")
                              message(cond); cat("\n")
                              return(c(1,1))
                            })
          if (length(model)==2) {
            next
          }

        } else if (model_choice=="gic.PDMR") {
          index=3
          if (run==1) {
            cat("PDMR method\n")
          }
          model <- tryCatch(DMRnet(insurance.train.10percent.x, insurance.train.10percent.y, algorithm="PDMR", family="gaussian"),
                            error=function(cond) {
                              message("Numerical instability in regular model creation detected. Will skip this set. Original error:")
                              message(cond); cat("\n")
                              return(c(1,1))
                            })
          if (length(model)==2) {
            next
          }

          gic <- gic.DMR(model, c = 2)   #we are using existing gic calculation which is compatible with PDMR models

        }  else  if (model_choice=="cv+sd.PDMR") {
          index=2
          if (run==1) {
            cat("PDMR with cv+sd\n")
          }
          model <- tryCatch(cv.DMRnet(insurance.train.10percent.x, insurance.train.10percent.y, family="gaussian", indexation.mode = "dimension", algorithm="PDMR"),
                            error=function(cond) {
                              message("Numerical instability in model creation in CV (cv.DMRnet) detected. Will skip this set. Original error:")
                              message(cond); cat("\n")
                              return(c(1,1))
                            })
          if (length(model)==2) {
            next
          }

        } else
          stop("Uknown method")

        if (model_choice=="gic.DMRnet" | model_choice=="gic.PDMR") {
          #cat(model_choice, "pred\n")
          prediction<- predict(model, newx=insurance.test.10percent.x, df = gic$df.min, type="class", unknown.factor.levels="NA")
          prediction1<- predict(gic, newx=insurance.test.10percent.x, type="class", unknown.factor.levels="NA")
          df<-gic$df.min
          cat("diff: ", sum(prediction1[!is.na(prediction1)]-prediction[!is.na(prediction)]), "\n")
        } else  if (model_choice=="cvg.DMRnet") {
          #cat(model_choice, "pred\n")
          prediction<- predict(model, newx=insurance.test.10percent.x, type="class", unknown.factor.levels="NA")
          df<-model$df.min
        } else  if (model_choice=="cv+sd.PDMR") {
          #cat(model_choice, "pred\n")
          prediction<- predict(model, newx=insurance.test.10percent.x, type="class", md="df.1se", unknown.factor.levels="NA")  #the very first test case of cv+sd.PDMR fails because of unknown factor levels, if not set to "NA"
          df<-model$df.1se
        } else
          stop("Uknown method")





        #prediction[is.na(prediction)] <- 0
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

postscript("result_insurance.ps", horizontal=TRUE,onefile=FALSE)
par(mfrow=c(2,4))

for (model_index in 1:4) {
  model_choice <- c("cv+sd.PDMR","gic.PDMR",  "cvg.DMRnet", "gic.DMRnet") [model_index]
  index <- model_index + 1 # a number of the corresponding expected result in a result filename
  vioplot(list(actual=results_mes[,model_index], expected=insurance.expected.errors[[index]][1:40]), ylab="Error", main=model_choice)
  vioplot(list(actual=results_dfs[,model_index], expected=insurance.expected.df[[index]][1:40]), ylab="Model Dimension", main=model_choice)
}


graphics.off()
