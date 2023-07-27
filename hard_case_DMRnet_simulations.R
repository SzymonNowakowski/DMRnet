library(devtools)
load_all()

#this file is the test case for issue#39 https://github.com/SzymonNowakowski/DMRnet/issues/39


all <- function(family, indexation.mode, algorithm) {

  print(">>>>> FAMILY:")
  print(family)
  print(">>>>> INDEXATION:")
  print(indexation.mode)
  print(">>>>> ALGO:")
  print(algorithm)


  #the below RData files automatically change .Random.seed because they contains this variable, but one cannot see it in the Environment in RStudio
  #to see it try to execute:
  # load("data/DMRnet_simulations/5_1_1_data.RData", envir <- new.env())
  # print(envir$.Random.seed)

    ##### THE BELOW 5 TEST CASES WERE ORIGINALLY IN MD-INDEXED CV, BUT HERE WE CHECK ALL COMBINATIONS ANYWAY

  print("5_1_1 seed+0, DMRnet")
  load("data/DMRnet_simulations/5_1_1_data.RData")
  if (family == "binomial") y <- factor(y>mean(y))
  model <- cv.DMRnet(XX, y, indexation.mode=indexation.mode, family=family, algorithm=algorithm, nlambda=20)  #in 0.3.2.9001 Error in 1:ncol(SS) : argument of length 0
                                                          # Called from: lapply(1:ncol(SS), function(i) DMRnet4lm_help(SS[, i], X, y,
                                                          #   fl, clust.method, lam))
  print("5_1_1 seed+1, DMRnet")
  load("data/DMRnet_simulations/5_1_1_plus_1_data.RData")
  if (family == "binomial") y <- factor(y>mean(y))
  model <- cv.DMRnet(XX, y, indexation.mode=indexation.mode, family=family, algorithm=algorithm, nlambda=20)  #in 0.3.2.9001 Error in SS[, i] : subscript out of bounds
                                                          # Called from: DMRnet4lm_help(SS[, i], X, y, fl, clust.method, lam)

  print("6_1_1 seed+0, DMRnet")
  load("data/DMRnet_simulations/6_1_1_data.RData")
  if (family == "binomial") y <- factor(y>mean(y))
  model <- cv.DMRnet(XX, y, indexation.mode=indexation.mode, family=family, algorithm=algorithm, nlambda=20)  #in 0.3.2.9001 Error in SS[, i] : subscript out of bounds
                                                          # Called from: DMRnet4lm_help(SS[, i], X, y, fl, clust.method, lam)

  print("6_1_1 seed+1, DMRnet")
  load("data/DMRnet_simulations/6_1_1_plus_1_data.RData")
  if (family == "binomial") y <- factor(y>mean(y))
  model <- cv.DMRnet(XX, y, indexation.mode=indexation.mode, family=family, algorithm=algorithm, nlambda=20)  #in 0.3.2.9001 Error in SS[, i] : subscript out of bounds
                                                          # Called from: DMRnet4lm_help(SS[, i], X, y, fl, clust.method, lam)


  print("6_2_1 seed+0, DMRnet (the fastest error happening in the full model)")
  load("data/DMRnet_simulations/6_2_1_data.RData")
  if (family == "binomial") y <- factor(y>mean(y))
  model <- cv.DMRnet(XX, y, indexation.mode=indexation.mode, family=family, algorithm=algorithm, nlambda=20)  #in 0.3.2.9001 Error in SS[, i] : subscript out of bounds
                                                          # Called from: DMRnet4lm_help(SS[, i], X, y, fl, clust.method, lam)
       # in correcting that case, I had to make a correction into cv_MD_indexed because foldmin was == 1
       # accordingly the test of plotting & prediction should be in order
  plot(model)
  #an error is expected here, but it should be a predefined DMRnet error, not some other failure
  tryCatch(predict(model, newx=head(XX), md="df.1se"), error=function(cond) {cat("Expected error: \n"); message(cond); cat("\n"); return(0)})
  tryCatch(coef(model, md="df.1se"), error=function(cond) {cat("Expected error: \n"); message(cond); cat("\n"); return(0)})


    ##### THE BELOW 3 TEST CASES WERE ORIGINALLY IN GIC-INDEXED CV, BUT HERE WE CHECK ALL COMBINATIONS ANYWAY


  print("6_2_1 seed+1, DMRnet")
  load("data/DMRnet_simulations/6_2_1_plus_1_data.RData")
  if (family == "binomial") y <- factor(y>mean(y))
  model <- cv.DMRnet(XX, y, indexation.mode=indexation.mode, family=family, algorithm=algorithm, nlambda=20)  #in 0.3.2.9001 Error in 1:ncol(SS) : argument of length 0
                             # Called from: lapply(1:ncol(SS), function(i) DMRnet4lm_help(SS[, i], X, y,
                             #   fl, clust.method, lam))

  print("6_3_1 seed+0, DMRnet")
  load("data/DMRnet_simulations/6_3_1_data.RData")
  if (family == "binomial") y <- factor(y>mean(y))
  model <- cv.DMRnet(XX, y, indexation.mode=indexation.mode, family=family, algorithm=algorithm, nlambda=20)  # in 0.3.2.9001 Error in 1:ncol(SS) : argument of length 0
                             # Called from: lapply(1:ncol(SS), function(i) DMRnet4lm_help(SS[, i], X, y,
                             # fl, clust.method, lam))

  print("4_2_1 seed+1, DMRnet")
  load("data/DMRnet_simulations/4_2_1_plus_1_data.RData")
  if (family == "binomial") y <- factor(y>mean(y))
  model <- cv.DMRnet(XX, y, indexation.mode=indexation.mode, family=family, algorithm=algorithm, nlambda=20)  # in 0.3.2.9001 Error in 1:ncol(SS) : argument of length 0
                             # Called from: lapply(1:ncol(SS), function(i) DMRnet4lm_help(SS[, i], X, y,
                             # fl, clust.method, lam))


}

all("gaussian", "dimension", "DMRnet")
all("gaussian", "GIC", "DMRnet")

all("gaussian", "dimension", "glamer")         # ?
all("gaussian", "GIC", "glamer")               # ?

all("gaussian", "dimension", "PDMR")         # ?
all("gaussian", "GIC", "PDMR")               # ?

all("binomial", "dimension", "DMRnet")         # it is not causing problems in 0.3.2.9001
all("binomial", "GIC", "DMRnet")               # it is not causing problems in 0.3.2.9001

all("binomial", "dimension", "glamer")         # ?
all("binomial", "GIC", "glamer")               # ?

all("binomial", "dimension", "PDMR")         # ?
all("binomial", "GIC", "PDMR")               # ?

#and finally a similarly caused error for SOSnet:
print("diagonal SOSnet with intercept")
DMRnet(diag(10), c(rep(0,5), rep(1,5)), interc=TRUE)   #it fails in 0.3.2.9001
                                                 #Error in bb[, i] : subscript out of bounds

print("diagonal SOSnet without intercept")
DMRnet(diag(10), c(rep(0,5), rep(1,5)), interc=FALSE)   #this doesn't cause problems in 0.3.2.9001


print("completed")
