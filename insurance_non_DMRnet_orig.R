library(digest)
library(grpreg)




run_list = c("scope", "scope", "lr", "cv.glmnet", "cv.MCP", "cv.grLasso")
runs<-100

part_filename_and_number = "insurance_orig"
print(part_filename_and_number)
set.seed(strtoi(substr(digest(part_filename_and_number, "md5", serialize = FALSE),1,7),16))
cat("seed set as md5 hash of the following string: ", part_filename_and_number,"\n\n")
cat("search for the result files postfixed with the same string: ", part_filename_and_number,"\n\n")



source("insurance_orig.R")
