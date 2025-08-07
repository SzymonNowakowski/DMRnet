library(devtools)
library(vioplot)
library(digest)
load_all()


#if no arguments, it generates data (it is a very long process, a couple of months at least on decent hardware) or load the precomputed results
#if arguments are provided (beta, snr, rho: in that order) it generates and saves the precomputed portion of data

##### Model matrix (data.frame) generation
generModel = function(n, p, rho){
  ## generation iid W[i,] ~ N(0,S), S=(1-a)*I_p + a*11'
  W <- matrix(rnorm(n*p),n,p)
  a <- 2*sin(pi/6*rho)
  b1 <- sqrt(1-a); b2 <- (sqrt(a*p+1-a)-b1)/p
  W <- t(apply(W,1,function(w) w <- b1*w+b2*sum(w) ) )
  #check:
  #U <- apply(W, 2, function(x) pnorm(x) ); print(round(100*cor(U)))
  XX <- lapply(data.frame(W), function(x) as.factor(ceiling(24*pnorm(x))) )
  as.data.frame(XX)
}









#####    M o d e l    p a r a m e t e r s
beta0 <- rep(0, 23)

# High 1 and 2
#betaHigh1.a <- c(rep(-2, 8),  rep(0, 8), rep(2, 8))
#betaHigh1.b <- c(rep(-2, 10), rep(0, 4), rep(2, 10))
betaHigh1.a <- c(rep(0, 8-1),  rep(2, 8), rep(4, 8))
betaHigh1.b <- c(rep(0, 10-1), rep(2, 4), rep(4, 10))
betaHigh1 <- c(-12, rep(betaHigh1.a, 3), rep(betaHigh1.b, 3), rep(beta0, 94) )

# High 3
betaHigh3.a <- c(rep(0, 8-1), rep(2, 8), rep(4, 8))
betaHigh3.b <- c(rep(0, 16-1), rep(5, 8))
betaHigh3 <- c(-14, rep(betaHigh3.a, 3), rep(betaHigh3.b, 3), rep(beta0, 94) )

# High 4
#betaHigh4.a <- c(rep(-2, 5),  rep(-1, 5), rep(0, 4), rep(1, 5), rep(2, 5))
betaHigh4.a <- c(rep(0, 5-1),  rep(1, 5), rep(2, 4), rep(3, 5), rep(4, 5))
betaHigh4 <- c(-10, rep(betaHigh4.a, 5), rep(beta0, 95))

# High 5 and 6
#betaHigh5.a <- c(rep(-2, 16), rep(3, 8))
betaHigh5.a <- c(rep(0, 16-1), rep(5, 8))
betaHigh5 <- c(-62, rep(betaHigh5.a, 25), rep(beta0, 75))

# High 7
#betaHigh7.a <- c(rep(-2, 4), rep(0, 12), rep(2, 8))
betaHigh7.a <- c(rep(0, 4-1), rep(2, 12), rep(4, 8))
betaHigh7 <- c(-20, rep(betaHigh7.a, 10), rep(beta0, 90) )

# High 8
#betaHigh8.a <- c(rep(-3, 6),  rep(-1, 6), rep(1, 6), rep(3, 6))
betaHigh8.a <- c(rep(0, 6-1),  rep(2, 6), rep(4, 6), rep(6, 6))
betaHigh8 <- c(-15, rep(betaHigh8.a, 5), rep(beta0, 95))

beta_list <- list(setHigh1 = betaHigh1, setHigh3 = betaHigh3, setHigh4 = betaHigh4, setHigh5 = betaHigh5, setHigh7 = betaHigh7, setHigh8 = betaHigh8)

#####    C h o i c e    o f   p a r a m e t e r s     o f     t h i s     run

n <- 500; p <- 100;

#####    P r e d i c t i o n    e r r o r    e s t i m a t i o n

gr <- rep(1:p, each=23)

args <- commandArgs(trailingOnly = TRUE)


if (length(args)==0) {
  pdf(file="result_high_dimensional_simulation.pdf", width=2100 / 25.4, height=2400 / 25.4, onefile=FALSE)   #units: inches, calculated from mm (2400x2100 in mm^2)
  par(mfrow=c(21,24))

}

if (length(args)>=1) {
  beta_choice <- as.numeric(args[1])
} else {
  beta_choice <- 1:6
}

if (length(args)>=2) {
  snr_choice <- 1.5^(-1:5)[as.numeric(args[2])]
} else {
  snr_choice <- 1.5^(-1:5)
}

if (length(args)>=3) {
  rho_choice <- c(0, 0.5)[as.numeric(args[3])]
} else {
  rho_choice <- c(0, 0.5)
}



for (actual_beta in beta_choice) {

  beta <- beta_list[[actual_beta]]       #6 choices
  denot <- names(beta_list)[actual_beta]

  for (rho in rho_choice) {
    for (snr in snr_choice) {
      theme <- theme <- paste(denot, snr, rho)
      seed <- strtoi(substr(digest(theme, "md5", serialize = FALSE),1,7),16)
      set.seed(seed)
      alg_list <- c("cvg.PDMR", "cv_sd.PDMR", "cvg.DMRnet")
      exp_list <- c("cvg.glamer", "cv_sd.glamer", "cvg.DMRnet")

      for (alg_no in 1:3) {
        alg <- alg_list[alg_no]

        cat("alg:", alg, "in", theme, "\n")

        filename <- paste(denot, "snr", as.character(round(snr,3)), "rho", as.character(rho), alg, sep="_")
        exp_filename <- paste(denot, "snr", as.character(round(snr,3)), "rho", as.character(rho), exp_list[alg_no], sep="_")

        filename_with_ext<-paste(filename,"csv",sep=".")
        exp_filename_with_ext<-paste(exp_filename,"csv",sep=".")

        expected_results<-read.table(paste("data_simulations/0.3.1/", exp_filename_with_ext, sep=""), header=TRUE, sep=",")
        rownames(expected_results)<-expected_results$X
        expected_results<-expected_results[-1]

        if (file.exists(filename_with_ext)) {
          OUT<-read.table(filename_with_ext, header=TRUE, sep=",")
          OUT<-OUT[,-1]
        } else {
          OUT <- simplify2array( lapply(1:200, function(i){
            cat(i,"\n")
            ##### 1. Model generation
            XX <- generModel(n, p, rho)
            X <- model.matrix(~., data=XX)
            mu <- X %*% beta
            signal2 <- mean( (mu - mean(mu))^2 )
            sigma <- sqrt(signal2) / snr
            y <- mu + rnorm(n, sd=sigma)

            ##### 2. Fitting methods

            MOD <- cv.DMRnet(XX, y, nlambda=100, nfolds=10, algorithm = ifelse(substr(alg, nchar(alg)-3, nchar(alg))=="PDMR", "PDMR", "DMRnet"), indexation.mode = (if(substr(alg,1,3)=="cvg") "GIC" else "dimension" ))
            md0 <- if(substr(alg,1,3)=="cvg") MOD$df.min else MOD$df.1se

            ##### 3. Prediction error estimation
            ERR <- replicate(100,{
              XX_test <- generModel(1000, p, rho)
              X_test <- model.matrix(~., data=XX_test)
              mu_test <- X_test %*% beta

              yy <- predict(MOD, XX_test)

              c( signal2_test = mean((mu_test - mean(mu_test))^2),
                 mse = mean((mu_test - yy)^2) )
            })
            c(signal2 = signal2, sigma = sigma, apply(ERR,1,mean), md0=md0)

          }))


          write.csv(OUT, file=filename_with_ext)
        }
        if (length(args)==0) {
          main_desc <- paste(denot, paste("snr", as.character(round(snr,3)), sep="="), paste("rho", as.character(rho), sep="="), sep=" ")

          vioplot(list(actual=as.numeric(OUT[4,]), expected=as.numeric(expected_results["mse",])), xlab = alg, ylab="Error", main=main_desc)
          if (sd(as.numeric(expected_results["md0",])) == 0 & sd(as.numeric(OUT[5,])) == 0)
            boxplot(list(actual=as.numeric(OUT[5,]), expected=as.numeric(expected_results["md0",])), xlab = alg, ylab="Model Dimension", main=main_desc)
          else
            vioplot(list(actual=as.numeric(OUT[5,]), expected=as.numeric(expected_results["md0",])), xlab = alg, ylab="Model Dimension", main=main_desc)
        }
      }
    }
  }
}


graphics.off()
