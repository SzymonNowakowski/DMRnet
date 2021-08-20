library(digest)
library(scriptName)
library(DMRnet)



run_list = c("cv.DMRnet", "gic.DMRnet")
runs<-50
sigma<-0.03

filename_and_dirs <- current_filename()
part_filename_and_number <- substr(filename_and_dirs, nchar(filename_and_dirs)-23, nchar(filename_and_dirs))
print(part_filename_and_number)
set.seed(strtoi(substr(digest(part_filename_and_number, "md5", serialize = FALSE),1,7),16))
cat("seed set as md5 hash of the following string: ", part_filename_and_number,"\n\n")
cat("search for the result files postfixed with the same string: ", part_filename_and_number,"\n\n")

source("insurance.R")
