library(devtools)
load_all()

stopper <- FALSE

###############################################################
###### Problems with constant intercept columns ###############
###############################################################



######### SCENARIO A SOSnet no intercept, matrix with intercept

cat("Sosnet\n")
X <- matrix(rnorm(50*15), ncol=15)
X[,1] <- 1.0    #constant
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)


mod <- DMRnet(X, y, interc=FALSE)   #SOSnet gets called - no factor columns
#pred <- predict(mod, newx=X, df=14)
MSEA <- sum((y - pred)^2)
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B SOSnet no intercept, matrix no interecept

X <- matrix(rnorm(50*15), ncol=15)
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)


mod <- DMRnet(X, y, interc=FALSE)   #SOSnet gets called - no factor columns
pred <- predict(mod, newx=X, df=15)
MSEB <- sum((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")

######### SCENARIO C SOSnet with intercept, matrix with intercept

X <- matrix(rnorm(50*15), ncol=15)
X[,1] <- 1.0    #constant
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)


mod <- DMRnet(X, y, interc=TRUE)   #SOSnet gets called - no factor columns
pred <- predict(mod, newx=X, df=15)
MSEC <- sum((y - pred)^2)
cat("Scenario C: ", MSEA, "\n")


######### SCENARIO D SOSnet with intercept, matrix no intercept

X <- matrix(rnorm(50*15), ncol=15)
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)


mod <- DMRnet(X, y, interc=TRUE)   #SOSnet gets called - no factor columns
#Warning message:
#  In cbind(be, c(mnk$coef, rep(0, ncol(X)))) :
#  number of rows of result is not a multiple of vector length (arg 2)
pred <- predict(mod, newx=X, df=16)
MSED <- sum((y - pred)^2)
cat("Scenario D: ", MSEB, "\n")



if (stopper)
  if ( MSEA > 1e-10 | MSEB > 1e-10 | MSEC > 1e-10 | MSED > 1e-10)
    stop("Incorrect computation with/without intercept in SOSnet")



###############################################################
###### Permutation of factor/numeric columns ##################
###############################################################



######### SCENARIO A DMR

cat("DMR\n")
X <- matrix(sample(1:2,750, replace=TRUE), ncol=15)
for (i in 1:13)
  X[,i] <- rnorm(50)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,14:15] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 14:15)
  X[,i] <- factor(X[,i])

mod <- DMR(X, y)
pred <- predict(mod, newx=X, df=3)
MSEA <- sum((y - pred)^2)
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B DMR

X <- matrix(sample(1:2,750, replace=TRUE), ncol=15)
for (i in 3:15)
  X[,i] <- rnorm(50)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,1:2] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 1:2)
  X[,i] <- factor(X[,i])

mod <- DMR(X, y)
pred <- predict(mod, newx=X, df=3)
MSEB <- sum((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")

if (stopper)
  if (MSEA / MSEB > 1000.0 | MSEA > 1e-10 | MSEB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))


######### SCENARIO A DMRnet

cat("DMRnet\n")
X <- matrix(sample(1:2,750, replace=TRUE), ncol=15)
for (i in 1:13)
  X[,i] <- rnorm(50)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,14:15] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 14:15)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y)
pred <- predict(mod, newx=X, df=3)
MSEA <- sum((y - pred)^2)
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B DMRnet

X <- matrix(sample(1:2,750, replace=TRUE), ncol=15)
for (i in 3:15)
  X[,i] <- rnorm(50)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,1:2] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 1:2)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y)
pred <- predict(mod, newx=X, df=3)
MSEB <- sum((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")


if (stopper)
  if (MSEA / MSEB > 1000.0 | MSEA > 1e-10 | MSEB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))



######### SCENARIO A GLAMER

cat("GLAMER\n")
X <- matrix(sample(1:2,750, replace=TRUE), ncol=15)
for (i in 1:13)
  X[,i] <- rnorm(50)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,14:15] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 14:15)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y, algorithm="glamer")
pred <- predict(mod, newx=X, df=3)
MSEA <- sum((y - pred)^2)
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B GLAMER

X <- matrix(sample(1:2,750, replace=TRUE), ncol=15)
for (i in 3:15)
  X[,i] <- rnorm(50)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,1:2] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 1:2)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y, algorithm="glamer")
pred <- predict(mod, newx=X, df=3)
MSEB <- sum((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")


if (stopper)
  if (MSEA / MSEB > 1000.0 | MSEA > 1e-10 | MSEB > 1e-10)
    stop(paste("Incorrect compuration for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))



cat("completed\n")
