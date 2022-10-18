library(devtools)
load_all()

#the below RData files automatically change .Random.seed because they contains this variable, but one cannot see it in the Environment in RStudio
#to see it try to execute:
# load("data/5_1_1_data.RData", envir <- new.env())
# print(envir$.Random.seed)

print("4_2_1 seed+1, DMRnet")
data("4_2_1_plus_1_data")
model <- cv.DMRnet(XX, y)  # in 0.3.2.9001 Error in 1:ncol(SS) : argument of length 0
                           # Called from: lapply(1:ncol(SS), function(i) DMRnet4lm_help(SS[, i], X, y,
                           # fl, clust.method, lam))



print("5_1_1 seed+0, DMRnet")
data("5_1_1_data")
model <- cv.DMRnet(XX, y, indexation.mode="dimension")  #in 0.3.2.9001 Error in 1:ncol(SS) : argument of length 0
                                                        # Called from: lapply(1:ncol(SS), function(i) DMRnet4lm_help(SS[, i], X, y,
                                                        #   fl, clust.method, lam))
print("5_1_1 seed+1, DMRnet")
data("5_1_1_plus_1_data")
model <- cv.DMRnet(XX, y, indexation.mode="dimension")  #in 0.3.2.9001 Error in SS[, i] : subscript out of bounds
                                                        # Called from: DMRnet4lm_help(SS[, i], X, y, fl, clust.method, lam)

print("6_1_1 seed+0, DMRnet")
data("6_1_1_data")
model <- cv.DMRnet(XX, y, indexation.mode="dimension")  #in 0.3.2.9001 Error in SS[, i] : subscript out of bounds
                                                        # Called from: DMRnet4lm_help(SS[, i], X, y, fl, clust.method, lam)

print("6_1_1 seed+1, DMRnet")
data("6_1_1_plus_1_data")
model <- cv.DMRnet(XX, y, indexation.mode="dimension")  #in 0.3.2.9001 Error in SS[, i] : subscript out of bounds
                                                        # Called from: DMRnet4lm_help(SS[, i], X, y, fl, clust.method, lam)


print("6_2_1 seed+0, DMRnet (the fastest error, probably in the 1st fold)")
data("6_2_1_data")
model <- cv.DMRnet(XX, y, indexation.mode="dimension")  #in 0.3.2.9001 Error in SS[, i] : subscript out of bounds
                                                        # Called from: DMRnet4lm_help(SS[, i], X, y, fl, clust.method, lam)


print("6_2_1 seed+1, DMRnet")
data("6_2_1_plus_1_data")
model <- cv.DMRnet(XX, y)  #in 0.3.2.9001 Error in 1:ncol(SS) : argument of length 0
                           # Called from: lapply(1:ncol(SS), function(i) DMRnet4lm_help(SS[, i], X, y,
                           #   fl, clust.method, lam))

print("6_3_1 seed+0, DMRnet")
data("6_3_1_data")
model <- cv.DMRnet(XX, y)  # in 0.3.2.9001 Error in 1:ncol(SS) : argument of length 0
                           # Called from: lapply(1:ncol(SS), function(i) DMRnet4lm_help(SS[, i], X, y,
                           # fl, clust.method, lam))
