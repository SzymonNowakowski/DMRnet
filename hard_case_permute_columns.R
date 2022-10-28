library(devtools)
load_all()

######### SCENARIO A

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


######### SCENARIO B

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


if (MSEA / MSEB > 1000.0 | MSEA > 1e-10 | MSEB > 1e-10)
  stop(paste("Incorrect compuation of permuted columns: MSE in scenario A =", MSEA, "while MSE in scenario B =", MSEB))

cat("completed\n")
