library(devtools)
load_all()
load("data/CV_airbnb/airbnb.train.RData")

cat("data loaded\n")

.Random.seed<-full_seed   #as per https://r-coder.com/set-seed-r/

cat("random seed set\n")

mod<-cv.DMRnet(data.train.percent.x, data.train.percent.y, algorithm="PDMR", indexation.mode = "dimension", nlambda=20)
#the problem in DMRnet v. 0.3.1 is that mod$df.min (which is the model dimension of the *best* model) is higher than
#ncol(mod$dmr.fit$beta) (which is equal to the largest model available)

cat("PDMR model computed for the previously identified hard case\n")

if (mod$df.min > ncol(mod$dmr.fit$beta))                              #in DMRnet v. 0.3.1: > mod$df.min
  stop("model dimension of the *best* model higher than the largest model available")    # [1] 55
                                                                                         # > ncol(mod$dmr.fit$beta)
                                                                                         # [1] 53

cat("Starting blind testing\n")

#now move to more extensive (and blind) tests
#first, define a gaussian function, which will *not* test model computation in each run
#(because it fails on LAPACK-related bugs, as is the case with insurance dataset)
#but rather will test the model prediction. It was what failed with CV in 0.3.1
gaussian <- function(allX, ally, factor_columns, model_choices, set_name, train_percent, runs) {

  errors<-list()

  set.seed(as.integer(10000*train_percent))

  for (model_choice in model_choices) {

    cat("now running: ", model_choice, "\n")

    MSEs<-MSEs.1se<-lengths<-rep(0,runs)

    run<-1

    while (run<=runs) {
      #cat("generating train/test sets\n")
      sample.percent <- sample(1:nrow(allX), train_percent*nrow(allX))
      data.train.percent.x <- allX[sample.percent, , drop=FALSE]
      data.train.percent.y <- ally[sample.percent,drop=FALSE]

      data.test.percent.x <- allX[-sample.percent, , drop=FALSE]
      data.test.percent.y <- ally[-sample.percent,drop=FALSE]


      #####RECOMPUTATION OF RELEVANT FACTORS in train set, to remove levels with no representative data (empty factors). Needed for random forest and glmnet
      for (i in factor_columns)
        data.train.percent.x[,i] <- factor(data.train.percent.x[,i])

      #remove data from test set with factors not present in train subsample as this causes predict() to fail
      for (i in factor_columns) {
        train.levels <- levels(data.train.percent.x[,i])
        data.test.percent.y<-data.test.percent.y[which(data.test.percent.x[,i] %in% train.levels)]
        data.test.percent.x<-data.test.percent.x[which(data.test.percent.x[,i] %in% train.levels),]
      }
      for (i in factor_columns)
        data.test.percent.x[,i] <- factor(data.test.percent.x[,i],levels = levels(data.train.percent.x[,i]) )   #recalculate factors now for new test


      #removing columns with only one value:
      singular_columns<-which(sapply(lapply(data.train.percent.x, unique), length)==1) #for continous columns length is 0
      if (length(singular_columns)>0) {
        data.test.percent.x <- data.test.percent.x[,-singular_columns]
        data.train.percent.x <- data.train.percent.x[,-singular_columns]
        #cat("removed", length(singular_columns), "columns due to singular values\n")
      }

      if (model_choice=="cvg.DMRnet" ) {
        cat("run", run, model_choice, "with cv indexed by gic\n")
        model.percent <- tryCatch(cv.DMRnet(data.train.percent.x, data.train.percent.y, nlambda=100, nfolds=10),
                                  error=function(cond) {
                                    message("Numerical instability in model creation in CV (cvg.DMRnet) detected. Will skip this 1-percent set. Original error:")
                                    message(cond); cat("\n")
                                    return(c(1,1))
                                  })

        if (length(model.percent)==2) {
          next
        }

      } else  if (model_choice=="cv.DMRnet" ) {
        cat("run", run, model_choice, "with cv indexed by model dimension\n")
        model.percent <- tryCatch(cv.DMRnet(data.train.percent.x, data.train.percent.y, nlambda=100, nfolds=10, indexation.mode="dimension"),
                                  error=function(cond) {
                                    message("Numerical instability in model creation in CV (cv.DMRnet) detected. Will skip this 1-percent set. Original error:")
                                    message(cond); cat("\n")
                                    return(c(1,1))
                                  })

        if (length(model.percent)==2) {
          next
        }

      } else  if (model_choice=="cvg.GLAMER") {
        cat("run", run, "GLAMER with cv indexed by GIC\n")
        model.percent <- tryCatch(cv.DMRnet(data.train.percent.x, data.train.percent.y, nlambda=100, nfolds=10, algorithm="glamer"),
                                  error=function(cond) {
                                    message("Numerical instability in model creation in CV (cvg.GLAMER) detected. Will skip this 1-percent set. Original error:")
                                    message(cond); cat("\n")
                                    return(c(1,1))
                                  })

        if (length(model.percent)==2) {
          next
        }


      } else  if (model_choice=="cv.GLAMER") {
        cat("run", run, "GLAMER with cv indexed by model dimension\n")
        model.percent <- tryCatch(cv.DMRnet(data.train.percent.x, data.train.percent.y, nlambda=100, nfolds=10, algorithm="glamer", indexation.mode="dimension"),
                                  error=function(cond) {
                                    message("Numerical instability in model creation in CV (cv.GLAMER) detected. Will skip this 1-percent set. Original error:")
                                    message(cond); cat("\n")
                                    return(c(1,1))
                                  })

        if (length(model.percent)==2) {
          next
        } else  if (model_choice=="cvg.PDMR") {
          cat("run", run, "PDMR with cv indexed by GIC\n")
          model.percent <- tryCatch(cv.DMRnet(data.train.percent.x, data.train.percent.y, nlambda=100, nfolds=10, algorithm="PDMR"),
                                    error=function(cond) {
                                      message("Numerical instability in model creation in CV (cvg.PDMR) detected. Will skip this 1-percent set. Original error:")
                                      message(cond); cat("\n")
                                      return(c(1,1))
                                    })

          if (length(model.percent)==2) {
            next
          }


        } else  if (model_choice=="cv.PDMR") {
          cat("run", run, "PDMR with cv indexed by model dimension\n")
          model.percent <- tryCatch(cv.DMRnet(data.train.percent.x, data.train.percent.y, nlambda=100, nfolds=10, algorithm="PDMR", indexation.mode="dimension"),
                                    error=function(cond) {
                                      message("Numerical instability in model creation in CV (cv.PDMR) detected. Will skip this 1-percent set. Original error:")
                                      message(cond); cat("\n")
                                      return(c(1,1))
                                    })

          if (length(model.percent)==2) {
            next
          }
      }

      prediction<- predict(model.percent, newx=data.test.percent.x, md="df.min")
      prediction.1se<- predict(model.percent, newx=data.test.percent.x, md="df.1se")

      MSEs[run]<-mean((prediction - data.test.percent.y)^2)
      MSEs.1se[run]<-mean((prediction.1se - data.test.percent.y)^2)

      if (is.na(median(MSEs[MSEs>0])))
        stop("Median for df.min model error is NA")
      if (is.na(median(MSEs.1se[MSEs.1se>0])))
        stop("Median for df.1se model error is NA")

      run<-run+1
    }

    cat(">>overall median df.min model error = ", median(MSEs), "\n")
    cat(">>overall median df.1se model error = ", median(MSEs.1se), "\n")

    ################## THE ACTUAL ERROR DETECTION HAPPENS HERE!############################
    if (is.na(median(MSEs)))
      stop("Median for df.min model error is NA")
    if (is.na(median(MSEs.1se)))
      stop("Median for df.1se model error is NA")

    errors[[model_choice]]<-MSEs
    errors[[paste(model_choice, ".1se", sep="")]]<-MSEs.1se

  }

  write.csv(errors, paste(set_name, train_percent, "errors.csv", sep="_"))   #the intention is too keep those results for future versions comparison

}

load("data_airbnb/airbnb.RData")   #TODO: the intention is to move that part into separate massive_airbnb.R tests later on
cat("data loaded\n")

model_choices<-c( "cvg.GLAMER","cv.GLAMER", "cvg.PDMR","cv.PDMR","cvg.DMRnet", "cv.DMRnet")

cat("running test for ", model_choices, "\n")

cat(">>starting 0.04 percent 50 runs test\n")
gaussian(air_X, air_y, factor_columns=1:4, model_choices=model_choices, set_name="airbnb", train_percent=0.04, runs=50)

cat("completed\n")
