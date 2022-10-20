library(devtools)
load_all()

#the below RData files automatically change .Random.seed because they contains this variable, but one cannot see it in the Environment in RStudio
#to see it try to execute:
# load("data/DMRnet_simulations/5_1_1_data.RData", envir <- new.env())
# print(envir$.Random.seed)

print("4_2_1 seed+1, DMRnet")
load("data/DMRnet_simulations/4_2_1_plus_1_data.RData")
model <- cv.DMRnet(XX, y)  # in 0.3.2.9001 Error in 1:ncol(SS) : argument of length 0
                           # Called from: lapply(1:ncol(SS), function(i) DMRnet4lm_help(SS[, i], X, y,
                           # fl, clust.method, lam))



print("5_1_1 seed+0, DMRnet")
load("data/DMRnet_simulations/5_1_1_data.RData")
model <- cv.DMRnet(XX, y, indexation.mode="dimension")  #in 0.3.2.9001 Error in 1:ncol(SS) : argument of length 0
                                                        # Called from: lapply(1:ncol(SS), function(i) DMRnet4lm_help(SS[, i], X, y,
                                                        #   fl, clust.method, lam))
print("5_1_1 seed+1, DMRnet")
load("data/DMRnet_simulations/5_1_1_plus_1_data.RData")
model <- cv.DMRnet(XX, y, indexation.mode="dimension")  #in 0.3.2.9001 Error in SS[, i] : subscript out of bounds
                                                        # Called from: DMRnet4lm_help(SS[, i], X, y, fl, clust.method, lam)

print("6_1_1 seed+0, DMRnet")
load("data/DMRnet_simulations/6_1_1_data.RData")
model <- cv.DMRnet(XX, y, indexation.mode="dimension")  #in 0.3.2.9001 Error in SS[, i] : subscript out of bounds
                                                        # Called from: DMRnet4lm_help(SS[, i], X, y, fl, clust.method, lam)

print("6_1_1 seed+1, DMRnet")
load("data/DMRnet_simulations/6_1_1_plus_1_data.RData")
model <- cv.DMRnet(XX, y, indexation.mode="dimension")  #in 0.3.2.9001 Error in SS[, i] : subscript out of bounds
                                                        # Called from: DMRnet4lm_help(SS[, i], X, y, fl, clust.method, lam)


print("6_2_1 seed+0, dm-indexed, DMRnet (the fastest error, probably in the 1st fold)")
load("data/DMRnet_simulations/6_2_1_data.RData")
model <- cv.DMRnet(XX, y, indexation.mode="dimension")  #in 0.3.2.9001 Error in SS[, i] : subscript out of bounds
                                                        # Called from: DMRnet4lm_help(SS[, i], X, y, fl, clust.method, lam)
     # in correcting that case, I had to make a correction into cv_MD_indexed because foldmin was == 1
     # accordingly the test of plotting should be in order
plot(model)
     # and of GIC-indexation:
print("6_2_1 seed+0, gic-indexed")
load("data/DMRnet_simulations/6_2_1_data.RData")
model <- cv.DMRnet(XX, y)                                #in 0.3.2.9001 Error in SS[, i] : subscript out of bounds
                                                         # Called from: DMRnet4lm_help(SS[, i], X, y, fl, clust.method, lam)

     # and of plotting:
plot(model)

print("6_2_1 seed+1, DMRnet")
load("data/DMRnet_simulations/6_2_1_plus_1_data.RData")
model <- cv.DMRnet(XX, y)  #in 0.3.2.9001 Error in 1:ncol(SS) : argument of length 0
                           # Called from: lapply(1:ncol(SS), function(i) DMRnet4lm_help(SS[, i], X, y,
                           #   fl, clust.method, lam))

print("6_3_1 seed+0, DMRnet")
load("data/DMRnet_simulations/6_3_1_data.RData")
model <- cv.DMRnet(XX, y)  # in 0.3.2.9001 Error in 1:ncol(SS) : argument of length 0
                           # Called from: lapply(1:ncol(SS), function(i) DMRnet4lm_help(SS[, i], X, y,
                           # fl, clust.method, lam))
