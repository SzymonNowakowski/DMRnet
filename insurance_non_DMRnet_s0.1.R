library(digest)




run_list = c("scope", "scope", "lr", "cv.glmnet")
runs<-1000
sigma<-0.1

part_filename_and_number = "insurance_sigma=0.1"
print(part_filename_and_number)
set.seed(strtoi(substr(digest(part_filename_and_number, "md5", serialize = FALSE),1,7),16))
cat("seed set as md5 hash of the following string: ", part_filename_and_number,"\n\n")
cat("search for the result files postfixed with the same string: ", part_filename_and_number,"\n\n")



source("insurance.R")
