library(devtools)
load_all()
data(airbnb.train)

cat("data loaded\n")

.Random.seed<-full_seed   #as per https://r-coder.com/set-seed-r/

cat("random seed set\n")

mod<-cv.DMRnet(data.train.percent.x, data.train.percent.y, algorithm="glamer", indexation.mode = "dimension")
#the problem in DMRnet v. 0.3.1 is that mod$df.min (which is the model dimension of the *best* model) is higher than
#ncol(mod$dmr.fit$beta) (which is equal to the largest model available)

cat("glamer model computed\n")

if (mod$df.min > ncol(mod$dmr.fit$beta))                              #in DMRnet v. 0.3.1: > mod$df.min
  stop("model dimension of the *best* model higher than the largest model available")    # [1] 55
                                                                                         # > ncol(mod$dmr.fit$beta)
                                                                                         # [1] 53

cat("completed\n")
