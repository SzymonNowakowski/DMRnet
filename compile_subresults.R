

compile <- function(directory, DMR_prefix, DMR_suffix, GLAF_prefix, all_prefix, all_errors, all_sizes, length = 100) {


  letters = c('1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k')

  errors<-read.csv(paste(directory, all_prefix, all_errors, "csv", sep=""), header=TRUE)
  sizes<-read.csv(paste(directory, all_prefix, all_sizes, "csv", sep=""), header=TRUE)

  DMRerrors<-matrix(nrow=0, ncol=2)
  for (letter in letters) {
    suberrors<-read.csv(paste(directory, DMR_prefix, letter, DMR_suffix, all_errors, "csv", sep=""), header=TRUE)
    DMRerrors <- rbind(DMRerrors, suberrors[,2:3])
  }

  DMRsizes<-matrix(nrow=0, ncol=2)
  for (letter in letters) {
    subsizes<-read.csv(paste(directory, DMR_prefix, letter, DMR_suffix, all_sizes, "csv", sep=""), header=TRUE)
    DMRsizes <- rbind(DMRsizes, subsizes[,2:3])
  }


  GLAFerrors<-matrix(nrow=0, ncol=2)
  for (letter in letters) {
    suberrors<-read.csv(paste(directory, GLAF_prefix, letter, DMR_suffix, all_errors, "csv", sep=""), header=TRUE)
    GLAFerrors <- rbind(GLAFerrors, suberrors[,2:3])
  }

  GLAFsizes<-matrix(nrow=0, ncol=2)
  for (letter in letters) {
    subsizes<-read.csv(paste(directory, GLAF_prefix, letter, DMR_suffix, all_sizes, "csv", sep=""), header=TRUE)
    GLAFsizes <- rbind(GLAFsizes, subsizes[,2:3])
  }
  errors_total <- cbind( GLAFerrors, DMRerrors, errors[1:length,2:8])
  sizes_total <- cbind( GLAFsizes, DMRsizes, sizes[1:length,2:8])

  write.csv(errors_total, paste(directory, all_prefix, "total_", all_errors, "csv", sep=""))
  write.csv(sizes_total, paste(directory, all_prefix, "total_", all_sizes, "csv", sep=""))

}

directory <- "results/"
DMR_prefix <- "insurance_DMRnet_orig"
DMR_suffix <- ".R_"
GLAF_prefix = "insurance_glamer_orig"
all_prefix <- "insurance_orig_"
all_errors <- "errors."
all_sizes <- "model_sizes."

compile(directory, DMR_prefix, DMR_suffix, GLAF_prefix, all_prefix, all_errors, all_sizes)
