library(devtools)
load_all()

stopper <- TRUE


set.seed(0)   #2 generates an error

##########################################################################
###### Problems with constant intercept columns - binomial ###############
##########################################################################





######### SCENARIO A SOSnet binomial no intercept, matrix no interecept

cat("Sosnet binomial\n")
X <- matrix(rnorm(500*15), ncol=15)
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)
y <- factor(y>mean(y))

mod <- DMRnet(X, y, interc=FALSE, family="binomial")   #SOSnet gets called - no factor columns
pred <- predict(mod, newx=X, df=15, type="class")
errorA <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario A: ", errorA, "\n")

######### SCENARIO B SOSnet binomial with intercept, matrix no intercept

X <- matrix(rnorm(500*15), ncol=15)
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)
y <- factor(y>mean(y))

mod <- DMRnet(X, y, interc=TRUE, family="binomial")   #SOSnet gets called - no factor columns
pred <- predict(mod, newx=X, df=16, type="class")
errorB <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario B: ", errorB, "\n")

if (stopper)
  if ( errorA > 1e-2 | errorB > 1e-2)
    stop("Incorrect computation with/without intercept in SOSnet")



#########################################################################################
###### Permutation of factor/numeric columns - binomial #################################
#########################################################################################



######### SCENARIO A DMR binomial

cat("DMR binomial\n")
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


######### SCENARIO B DMR binomial

X <- X[,c(14, 15, 1:13)]

mod <- DMR(X, y, family="binomial")
pred <- predict(mod, newx=X, df=7, type="class")
errorB <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario B: ", errorB, "\n")

if (stopper)
  if (errorA!=errorB | errorA > 1e-10 | errorB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: error in scenario A =", errorA, "while error in scenario B =", errorB))


######### SCENARIO A DMRnet binomial

cat("DMRnet binomial\n")
X <- X[,c(3:15, 1, 2)]

mod <- DMRnet(X, y, family="binomial")
pred <- predict(mod, newx=X, df=7, type="class")
errorA <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario A: ", errorA, "\n")


######### SCENARIO B DMRnet binomial

X <- X[,c(14, 15, 1:13)]

mod <- DMRnet(X, y, family="binomial")
pred <- predict(mod, newx=X, df=7, type="class")
errorB <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario B: ", errorB, "\n")


if (stopper)
  if (errorA!=errorB | errorA > 1e-10 | errorB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: error in scenario A =", errorA, "while error in scenario B =", errorB))



######### SCENARIO A GLAMER binomial

cat("GLAMER binomial\n")
X <- X[,c(3:15, 1, 2)]

mod <- DMRnet(X, y, algorithm="glamer", family="binomial")
pred <- predict(mod, newx=X, df=7, type="class")                        #WHY df=2 for 2 x 2-factor columns???
errorA <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario A: ", errorA, "\n")


######### SCENARIO B GLAMER binomial

X <- X[,c(14, 15, 1:13)]

mod <- DMRnet(X, y, algorithm="glamer", family="binomial")
pred <- predict(mod, newx=X, df=7, type="class")                   #WHY df=2 for 2 x 2-factor columns???
errorB <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario B: ", errorB, "\n")


if (stopper)
  if (errorA!=errorB | errorA > 1e-10 | errorB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: error in scenario A =", errorA, "while error in scenario B =", errorB))


######### SCENARIO A PDMR binomial

cat("PDMR binomial\n")
X <- X[,c(3:15, 1, 2)]

mod <- DMRnet(X, y, algorithm="PDMR", family="binomial")
pred <- predict(mod, newx=X, df=7, type="class")                        #WHY df=2 for 2 x 2-factor columns???
errorA <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario A: ", errorA, "\n")


######### SCENARIO B PDMR binomial

X <- X[,c(14, 15, 1:13)]

mod <- DMRnet(X, y, algorithm="PDMR", family="binomial")
pred <- predict(mod, newx=X, df=7, type="class")                   #WHY df=2 for 2 x 2-factor columns???
errorB <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario B: ", errorB, "\n")


if (stopper)
  if (errorA!=errorB | errorA > 1e-10 | errorB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: error in scenario A =", errorA, "while error in scenario B =", errorB))

#########################################################################################
###### Permutation of factor/numeric columns - binomial with p>n ########################
#########################################################################################




######### SCENARIO A DMRnet binomial p>n

cat("DMRnet binomial p>n\n")
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


######### SCENARIO B DMRnet binomial p>n

X <- X[,c(14, 15, 1:13)]

mod <- DMRnet(X, y, family="binomial")
pred <- predict(mod, newx=X, df=4, type="class")
errorB <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario B: ", errorB, "\n")


if (stopper)
  if (errorA!=errorB | errorA > 1e-10 | errorB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: error in scenario A =", errorA, "while error in scenario B =", errorB))



######### SCENARIO A GLAMER binomial p>n

cat("GLAMER binomial p>n\n")
X <- X[,c(3:15, 1, 2)]

mod <- DMRnet(X, y, algorithm="glamer", family="binomial")
pred <- predict(mod, newx=X, df=4, type="class")                        #WHY df=2 for 2 x 2-factor columns???
errorA <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario A: ", errorA, "\n")


######### SCENARIO B GLAMER binomial p>n

X <- X[,c(14, 15, 1:13)]

mod <- DMRnet(X, y, algorithm="glamer", family="binomial")
pred <- predict(mod, newx=X, df=4, type="class")                   #WHY df=2 for 2 x 2-factor columns???
errorB <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario B: ", errorB, "\n")


if (stopper)
  if (errorA!=errorB | errorA > 1e-10 | errorB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: error in scenario A =", errorA, "while error in scenario B =", errorB))


######### SCENARIO A PDMR binomial p>n

cat("PDMR binomial p>n\n")
X <- X[,c(3:15, 1, 2)]

mod <- DMRnet(X, y, algorithm="PDMR", family="binomial")
pred <- predict(mod, newx=X, df=4, type="class")                        #WHY df=2 for 2 x 2-factor columns???
errorA <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario A: ", errorA, "\n")


######### SCENARIO B PDMR binomial p>n

X <- X[,c(14, 15, 1:13)]

mod <- DMRnet(X, y, algorithm="PDMR", family="binomial")
pred <- predict(mod, newx=X, df=4, type="class")                   #WHY df=2 for 2 x 2-factor columns???
errorB <- mean(as.integer(y) != (as.integer(pred)+1))
cat("Scenario B: ", errorB, "\n")


if (stopper)
  if (errorA!=errorB | errorA > 1e-10 | errorB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: error in scenario A =", errorA, "while error in scenario B =", errorB))




##########################################################################
###### Problems with constant intercept columns - gaussian ###############
##########################################################################




######### SCENARIO A SOSnet gaussian no intercept, matrix no interecept
cat("Sosnet\n")
X <- matrix(rnorm(500*15), ncol=15)
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)

mod <- DMRnet(X, y, interc=FALSE)   #SOSnet gets called - no factor columns
pred <- predict(mod, newx=X, df=15)
MSEA <- mean((y - pred)^2)
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B SOSnet gaussian with intercept, matrix no intercept

X <- matrix(rnorm(500*15), ncol=15)
b_vector <- sample(1:10,15, replace=TRUE)
y <- X %*% matrix(b_vector)

mod <- DMRnet(X, y, interc=TRUE)   #SOSnet gets called - no factor columns
pred <- predict(mod, newx=X, df=16)
MSEB <- mean((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")

if (stopper)
  if (MSEA > 1e-10 | MSEB > 1e-10)
    stop("Incorrect computation with/without intercept in SOSnet")



#########################################################################################
###### Permutation of factor/numeric columns - gaussian #################################
#########################################################################################



######### SCENARIO A DMR gaussian

cat("DMR gaussian\n")
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

######### SCENARIO B DMR gaussian

X <- X[,c(14, 15, 1:13)]

mod <- DMR(X, y)
pred <- predict(mod, newx=X, df=7)
MSEB <- mean((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")

if (stopper)
  if (MSEA != MSEB |MSEA > 1e-10 | MSEB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))


######### SCENARIO A DMRnet gaussian

cat("DMRnet gaussian\n")
X <- X[,c(3:15, 1, 2)]

mod <- DMRnet(X, y)
pred <- predict(mod, newx=X, df=7)
MSEA <- mean((y - pred)^2)
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B DMRnet gaussian

X <- X[,c(14, 15, 1:13)]

mod <- DMRnet(X, y)
pred <- predict(mod, newx=X, df=7)
MSEB <- mean((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")


if (stopper)
  if (MSEA != MSEB | MSEA > 1e-10 | MSEB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))



######### SCENARIO A GLAMER gaussian

cat("GLAMER gaussian\n")
X <- X[,c(3:15, 1, 2)]

mod <- DMRnet(X, y, algorithm="glamer")
pred <- predict(mod, newx=X, df=7)
MSEA <- mean((y - pred)^2)
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B GLAMER gaussian

X <- X[,c(14, 15, 1:13)]

mod <- DMRnet(X, y, algorithm="glamer")
pred <- predict(mod, newx=X, df=7)
MSEB <- mean((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")


if (stopper)
  if (MSEA != MSEB | MSEA > 1e-10 | MSEB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))

######### SCENARIO A PDMR gaussian

cat("PDMR gaussian\n")
X <- X[,c(3:15, 1, 2)]

mod <- DMRnet(X, y, algorithm="PDMR")
pred <- predict(mod, newx=X, df=7)
MSEA <- mean((y - pred)^2)
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B PDMR gaussian

X <- X[,c(14, 15, 1:13)]

mod <- DMRnet(X, y, algorithm="PDMR")
pred <- predict(mod, newx=X, df=7)
MSEB <- mean((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")


if (stopper)
  if (MSEA != MSEB | MSEA > 1e-10 | MSEB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))



#########################################################################################
###### Permutation of factor/numeric columns - gaussian with p>n ########################
#########################################################################################



######### SCENARIO A DMRnet gaussian p>n

cat("DMRnet gaussian p>n\n")
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


######### SCENARIO B DMRnet gaussian p>n

X <- X[,c(14, 15, 1:13)]

mod <- DMRnet(X, y)
pred <- predict(mod, newx=X, df=6)
MSEB <- mean((y - pred)^2)
cat("Scenario B: ", MSEB, "\n")


if (stopper)
  if (MSEA != MSEB | MSEA > 1e-10 | MSEB > 1e-10)
    stop(paste("Incorrect computation for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))



######### SCENARIO A GLAMER gaussian p>n

cat("GLAMER gaussian p>n\n")
X <- X[,c(3:15, 1, 2)]

mod <- DMRnet(X, y, algorithm="glamer")
pred <- predict(mod, newx=X, df=6)
MSEA <- mean((y - pred)^2)                #the MSE tends to be > 1.0 in GLAMER
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B GLAMER gaussian p>n

X <- X[,c(14, 15, 1:13)]

mod <- DMRnet(X, y, algorithm="glamer")
pred <- predict(mod, newx=X, df=6)
MSEB <- mean((y - pred)^2)                #the MSE tends to be > 1.0 in GLAMER
cat("Scenario B: ", MSEB, "\n")


if (stopper)
  if (abs(MSEA - MSEB) > 1e-10 )
    stop(paste("Incorrect computation for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))


######### SCENARIO A PDMR gaussian p>n

cat("PDMR gaussian p>n\n")
X <- X[,c(3:15, 1, 2)]

mod <- DMRnet(X, y, algorithm="PDMR")
pred <- predict(mod, newx=X, df=6)
MSEA <- mean((y - pred)^2)                #the MSE tends to be > 1.0 in PDMR
cat("Scenario A: ", MSEA, "\n")


######### SCENARIO B PDMR gaussian p>n

X <- X[,c(14, 15, 1:13)]

mod <- DMRnet(X, y, algorithm="PDMR")
pred <- predict(mod, newx=X, df=6)
MSEB <- mean((y - pred)^2)                #the MSE tends to be > 1.0 in PDMR
cat("Scenario B: ", MSEB, "\n")


if (stopper)
  if (abs(MSEA - MSEB) > 1e-10)
    stop(paste("Incorrect computation for permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))



cat("completed\n")
