library(devtools)
load_all()

stopper <- FALSE


set.seed(0)

##########################################################################
###### Problems with constant intercept columns - binomial ###############
##########################################################################


######### SCENARIO A SOSnet no intercept, matrix with intercept

cat("Sosnet\n")
X <- matrix(rnorm(500*15), ncol=15)
X[,1] <- 1.0    #constant
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)
y <- factor(y>mean(y))

mod <- DMRnet(X, y, interc=FALSE, family="binomial")   #SOSnet gets called - no factor columns
#pred <- predict(mod, newx=X, df=14, type="class")
pred<-rnorm(length(y))
errorA <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario A: ", errorA, "\n")


######### SCENARIO B SOSnet no intercept, matrix no interecept

X <- matrix(rnorm(500*15), ncol=15)
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)
y <- factor(y>mean(y))

mod <- DMRnet(X, y, interc=FALSE, family="binomial")   #SOSnet gets called - no factor columns
pred <- predict(mod, newx=X, df=15, type="class")
errorB <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario B: ", errorB, "\n")

######### SCENARIO C SOSnet with intercept, matrix with intercept

X <- matrix(rnorm(500*15), ncol=15)
X[,1] <- 1.0    #constant
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)
y <- factor(y>mean(y))

mod <- DMRnet(X, y, interc=TRUE, family="binomial")   #SOSnet gets called - no factor columns
#pred <- predict(mod, newx=X, df=15, type="class")
pred<-rnorm(length(y))
errorC <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario C: ", errorA, "\n")


######### SCENARIO D SOSnet with intercept, matrix no intercept

X <- matrix(rnorm(500*15), ncol=15)
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)
y <- factor(y>mean(y))

mod <- DMRnet(X, y, interc=TRUE, family="binomial")   #SOSnet gets called - no factor columns
pred <- predict(mod, newx=X, df=16, type="class")
errorD <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario D: ", errorB, "\n")

if (stopper)
  if ( errorA > 1e-10 | errorB > 1e-10 | errorC > 1e-10 | errorD > 1e-10)
    stop("Incorrect computation with/without intercept in SOSnet")



#########################################################################################
###### Permutation of factor/numeric columns - binomial #################################
#########################################################################################



######### SCENARIO A DMR

cat("DMR\n")
X <- matrix(sample(1:4,7500, replace=TRUE), ncol=15)
for (i in 1:13)
  X[,i] <- rnorm(500)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,14:15] %*% matrix(b_vector)
y <- factor(y>mean(y))

X<-data.frame(X)
for (i in 14:15)
  X[,i] <- factor(X[,i])

mod <- DMR(X, y, family="binomial")
pred <- predict(mod, newx=X, df=7, type="class")
errorA <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario A: ", errorA, "\n")


######### SCENARIO B DMR

X <- matrix(sample(1:4,7500, replace=TRUE), ncol=15)
for (i in 3:15)
  X[,i] <- rnorm(500)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,1:2] %*% matrix(b_vector)
y <- factor(y>mean(y))

X<-data.frame(X)
for (i in 1:2)
  X[,i] <- factor(X[,i])

mod <- DMR(X, y, family="binomial")
pred <- predict(mod, newx=X, df=7, type="class")
errorB <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario B: ", errorB, "\n")

if (stopper)
  if (errorA / errorB > 1000.0 | errorA > 1e-10 | errorB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: error in scenario A =", errorA, "while error in scenario B =", errorB))


######### SCENARIO A DMRnet

cat("DMRnet\n")
X <- matrix(sample(1:4,7500, replace=TRUE), ncol=15)
for (i in 1:13)
  X[,i] <- rnorm(500)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,14:15] %*% matrix(b_vector)
y <- factor(y>mean(y))

X<-data.frame(X)
for (i in 14:15)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y, family="binomial")
pred <- predict(mod, newx=X, df=7, type="class")
errorA <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario A: ", errorA, "\n")


######### SCENARIO B DMRnet

X <- matrix(sample(1:4,7500, replace=TRUE), ncol=15)
for (i in 3:15)
  X[,i] <- rnorm(500)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,1:2] %*% matrix(b_vector)
y <- factor(y>mean(y))

X<-data.frame(X)
for (i in 1:2)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y, family="binomial")
pred <- predict(mod, newx=X, df=7, type="class")
errorB <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario B: ", errorB, "\n")


if (stopper)
  if (errorA / errorB > 1000.0 | errorA > 1e-10 | errorB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: error in scenario A =", errorA, "while error in scenario B =", errorB))



######### SCENARIO A GLAMER

cat("GLAMER\n")
X <- matrix(sample(1:4,7500, replace=TRUE), ncol=15)
for (i in 1:13)
  X[,i] <- rnorm(500)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,14:15] %*% matrix(b_vector)
y <- factor(y>mean(y))

X<-data.frame(X)
for (i in 14:15)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y, algorithm="glamer", family="binomial")
pred <- predict(mod, newx=X, df=7, type="class")                        #WHY df=2 for 2 x 2-factor columns???
errorA <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario A: ", errorA, "\n")


######### SCENARIO B GLAMER

X <- matrix(sample(1:4,7500, replace=TRUE), ncol=15)
for (i in 3:15)
  X[,i] <- rnorm(500)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,1:2] %*% matrix(b_vector)
y <- factor(y>mean(y))

X<-data.frame(X)
for (i in 1:2)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y, algorithm="glamer", family="binomial")
pred <- predict(mod, newx=X, df=7, type="class")                   #WHY df=2 for 2 x 2-factor columns???
errorB <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario B: ", errorB, "\n")


if (stopper)
  if (errorA / errorB > 1000.0 | errorA > 1e-10 | errorB > 1e-10)
    stop(paste("Incorrect compuration for permuted columns: error in scenario A =", errorA, "while error in scenario B =", errorB))

#########################################################################################
###### Permutation of factor/numeric columns - binomial with p>n ########################
#########################################################################################




######### SCENARIO A DMRnet

cat("DMRnet p>n\n")
X <- matrix(sample(1:4,150, replace=TRUE), ncol=15)
for (i in 1:13)
  X[,i] <- rnorm(10)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,14:15] %*% matrix(b_vector)
y <- factor(y>mean(y))

X<-data.frame(X)
for (i in 14:15)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y, family="binomial")
pred <- predict(mod, newx=X, df=4, type="class")
errorA <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario A: ", errorA, "\n")


######### SCENARIO B DMRnet

X <- matrix(sample(1:4,150, replace=TRUE), ncol=15)
for (i in 3:15)
  X[,i] <- rnorm(10)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,1:2] %*% matrix(b_vector)
y <- factor(y>mean(y))

X<-data.frame(X)
for (i in 1:2)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y, family="binomial")
pred <- predict(mod, newx=X, df=4, type="class")
errorB <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario B: ", errorB, "\n")


if (stopper)
  if (errorA / errorB > 1000.0 | errorA > 1e-10 | errorB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: error in scenario A =", errorA, "while error in scenario B =", errorB))



######### SCENARIO A GLAMER

cat("GLAMER p>n\n")
X <- matrix(sample(1:4,150, replace=TRUE), ncol=15)
for (i in 1:13)
  X[,i] <- rnorm(10)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,14:15] %*% matrix(b_vector)
y <- factor(y>mean(y))

X<-data.frame(X)
for (i in 14:15)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y, algorithm="glamer", family="binomial")
pred <- predict(mod, newx=X, df=4, type="class")                        #WHY df=2 for 2 x 2-factor columns???
errorA <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario A: ", errorA, "\n")


######### SCENARIO B GLAMER

X <- matrix(sample(1:4,150, replace=TRUE), ncol=15)
for (i in 3:15)
  X[,i] <- rnorm(10)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,1:2] %*% matrix(b_vector)
y <- factor(y>mean(y))

X<-data.frame(X)
for (i in 1:2)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y, algorithm="glamer", family="binomial")
pred <- predict(mod, newx=X, df=4, type="class")                   #WHY df=2 for 2 x 2-factor columns???
errorB <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario B: ", errorB, "\n")


if (stopper)
  if (errorA / errorB > 1000.0 | errorA > 1e-10 | errorB > 1e-10)
    stop(paste("Incorrect compuration for permuted columns: error in scenario A =", errorA, "while error in scenario B =", errorB))





##########################################################################
###### Problems with constant intercept columns - gaussian ###############
##########################################################################


######### SCENARIO A SOSnet no intercept, matrix with intercept

cat("Sosnet\n")
X <- matrix(rnorm(500*15), ncol=15)
X[,1] <- 1.0    #constant
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)

mod <- DMRnet(X, y, interc=FALSE)   #SOSnet gets called - no factor columns
#pred <- predict(mod, newx=X, df=14)
pred<-rnorm(length(y))
MSEA <- mean((y - pred)^2)
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B SOSnet no intercept, matrix no interecept

X <- matrix(rnorm(500*15), ncol=15)
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)

mod <- DMRnet(X, y, interc=FALSE)   #SOSnet gets called - no factor columns
pred <- predict(mod, newx=X, df=15)
MSEB <- mean((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")

######### SCENARIO C SOSnet with intercept, matrix with intercept

X <- matrix(rnorm(500*15), ncol=15)
X[,1] <- 1.0    #constant
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)

mod <- DMRnet(X, y, interc=TRUE)   #SOSnet gets called - no factor columns
#pred <- predict(mod, newx=X, df=15)
pred<-rnorm(length(y))
MSEC <- mean((y - pred)^2)
cat("Scenario C: ", MSEA, "\n")


######### SCENARIO D SOSnet with intercept, matrix no intercept

X <- matrix(rnorm(500*15), ncol=15)
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)

mod <- DMRnet(X, y, interc=TRUE)   #SOSnet gets called - no factor columns
pred <- predict(mod, newx=X, df=16)
MSED <- mean((y - pred)^2)
cat("Scenario D: ", MSEB, "\n")

if (stopper)
  if ( MSEA > 1e-10 | MSEB > 1e-10 | MSEC > 1e-10 | MSED > 1e-10)
    stop("Incorrect computation with/without intercept in SOSnet")



#########################################################################################
###### Permutation of factor/numeric columns - gaussian #################################
#########################################################################################



######### SCENARIO A DMR

cat("DMR\n")
X <- matrix(sample(1:4,7500, replace=TRUE), ncol=15)
for (i in 1:13)
  X[,i] <- rnorm(500)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,14:15] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 14:15)
  X[,i] <- factor(X[,i])

mod <- DMR(X, y)
pred <- predict(mod, newx=X, df=7)
MSEA <- mean((y - pred)^2)
cat("Scenario A: ", MSEA, "\n")

######### SCENARIO B DMR

X <- matrix(sample(1:4,7500, replace=TRUE), ncol=15)
for (i in 3:15)
  X[,i] <- rnorm(500)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,1:2] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 1:2)
  X[,i] <- factor(X[,i])

mod <- DMR(X, y)
pred <- predict(mod, newx=X, df=7)
MSEB <- mean((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")

if (stopper)
  if (MSEA / MSEB > 1000.0 | MSEA > 1e-10 | MSEB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))


######### SCENARIO A DMRnet

cat("DMRnet\n")
X <- matrix(sample(1:4,7500, replace=TRUE), ncol=15)
for (i in 1:13)
  X[,i] <- rnorm(500)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,14:15] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 14:15)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y)
pred <- predict(mod, newx=X, df=7)
MSEA <- mean((y - pred)^2)
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B DMRnet

X <- matrix(sample(1:4,7500, replace=TRUE), ncol=15)
for (i in 3:15)
  X[,i] <- rnorm(500)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,1:2] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 1:2)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y)
pred <- predict(mod, newx=X, df=7)
MSEB <- mean((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")


if (stopper)
  if (MSEA / MSEB > 1000.0 | MSEA > 1e-10 | MSEB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))



######### SCENARIO A GLAMER

cat("GLAMER\n")
X <- matrix(sample(1:4,7500, replace=TRUE), ncol=15)
for (i in 1:13)
  X[,i] <- rnorm(500)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,14:15] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 14:15)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y, algorithm="glamer")
pred <- predict(mod, newx=X, df=7)
MSEA <- mean((y - pred)^2)
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B GLAMER

X <- matrix(sample(1:4,7500, replace=TRUE), ncol=15)
for (i in 3:15)
  X[,i] <- rnorm(500)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,1:2] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 1:2)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y, algorithm="glamer")
pred <- predict(mod, newx=X, df=7)
MSEB <- mean((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")


if (stopper)
  if (MSEA / MSEB > 1000.0 | MSEA > 1e-10 | MSEB > 1e-10)
    stop(paste("Incorrect compuration for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))



#########################################################################################
###### Permutation of factor/numeric columns - gaussian with p>n ########################
#########################################################################################



######### SCENARIO A DMRnet

cat("DMRnet p>n\n")
X <- matrix(sample(1:4,150, replace=TRUE), ncol=15)
for (i in 1:13)
  X[,i] <- rnorm(10)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,14:15] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 14:15)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y)
pred <- predict(mod, newx=X, df=6)
MSEA <- mean((y - pred)^2)
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B DMRnet

X <- matrix(sample(1:4,150, replace=TRUE), ncol=15)
for (i in 3:15)
  X[,i] <- rnorm(10)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,1:2] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 1:2)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y)
pred <- predict(mod, newx=X, df=6)
MSEB <- mean((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")


if (stopper)
  if (MSEA / MSEB > 1000.0 | MSEA > 1e-10 | MSEB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))



######### SCENARIO A GLAMER

cat("GLAMER p>n\n")
X <- matrix(sample(1:4,150, replace=TRUE), ncol=15)
for (i in 1:13)
  X[,i] <- rnorm(10)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,14:15] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 14:15)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y, algorithm="glamer")
pred <- predict(mod, newx=X, df=6)
MSEA <- mean((y - pred)^2)
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B GLAMER

X <- matrix(sample(1:4,150, replace=TRUE), ncol=15)
for (i in 3:15)
  X[,i] <- rnorm(10)
b_vector <- sample(1:10,2, replace=TRUE)
y <- X[,1:2] %*% matrix(b_vector)

X<-data.frame(X)
for (i in 1:2)
  X[,i] <- factor(X[,i])

mod <- DMRnet(X, y, algorithm="glamer")
pred <- predict(mod, newx=X, df=6)
MSEB <- mean((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")


if (stopper)
  if (MSEA / MSEB > 1000.0 | MSEA > 1e-10 | MSEB > 1e-10)
    stop(paste("Incorrect compuration for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))



cat("completed\n")
