cv_helper<-function(Xtr, ytr, Xte, yte, real_n, agressive) {


  ####SzN remove from train and test columns causing data singularity
  #removing columns with only one level:
  singular_factors<-which(sapply(lapply(Xtr, levels), length)==1)  #for continous columns length is 0
  if (length(singular_factors)>0) {
    Xte <- Xte[,-singular_factors]
    Xtr <- Xtr[,-singular_factors]
    cat("removed", length(singular_factors), "columns due to singular factors in training set\n")
  }

  ###SzN remove from test the data with factors not present in training
  nn <- sapply(1:ncol(Xte), function(i) class(Xte[,i]))
  faki <- which(nn == "factor")
  n.factors <- length(faki)
  if (n.factors > 0)
    for (i in 1:n.factors) {
      Xtr[,faki[i]] <- factor(Xtr[,faki[i]])  #may be removed in next version of package because
      train.levels <- levels(Xtr[,faki[i]])   #the same info is in dmr$levels.listed[[i]]
      #but this has the advantage that is also compatible with the old package version
      yte<-yte[which(Xte[,faki[i]] %in% train.levels)]
      Xte<-Xte[which(Xte[,faki[i]] %in% train.levels),]
      Xte[,faki[i]]<-factor(Xte[,faki[i]], levels=train.levels)  #recalculate the factors in new test set, may be removed in next version of package because the same happens in predict(..)
    }
  real_n <- real_n + length(yte)
  #TODO: what to do if all test data is removed?


  if (agressive) {
    #preparation to detect data dependency
    Xtr.make <- makeX(Xtr)
    prev_pos <- 0
    first_identified_yet = FALSE
    for (i in 1:ncol(Xtr))
      if (i %in% faki) {  #removing columns from the last level (but not for the first original column) they are linearly dependant
        # cat(i, prev_pos, length(levels(insurance.train.10percent.x[,i])), "\n")
        if (first_identified_yet) {
          Xtr.make <- Xtr.make[,-(prev_pos+length(levels(Xtr[,i])))]
          prev_pos <- prev_pos+length(levels(Xtr[,i])) - 1
        } else {
          prev_pos <- prev_pos+length(levels(Xtr[,i]))
          first_identified_yet = TRUE
        }
      } else prev_pos<-prev_pos+1

    QR<- qr(Xtr.make)

    if (QR$rank < ncol(Xtr.make)) {  #singular
      reverse_lookup<-rep(0, ncol(Xtr.make))
      pos<-1
      first_identified_yet = FALSE
      for (i in 1:ncol(Xtr))
        if (i %in% faki) {
          if (first_identified_yet) {
            reverse_lookup[pos:(pos+length(levels(Xtr[,i]))-2)]<-i  #there are levels-1 columns corresponding to the first original column
            pos<-pos+length(levels(Xtr[,i]))-1
          } else {
            reverse_lookup[pos:(pos+length(levels(Xtr[,i]))-1)]<-i  #there are levels-1 columns corresponding to each original column other than the first
            pos<-pos+length(levels(Xtr[,i]))
            first_identified_yet = TRUE
          }
        } else {
          reverse_lookup[pos]<-i
          pos<-pos+1
        }

      #removal of columns for pivot positions larger than rank
      remove_us<-reverse_lookup[QR$pivot[(QR$rank+1):length(QR$pivot)]]
      Xtr <- Xtr[,-remove_us]
      Xte <- Xte[,-remove_us]
      cat("removed", length(unique(remove_us)), "columns due to data linear dependency\n")
    }
  } else
    cat("non-aggressive mode for singularity detection, only single-leveled factors handled\n")

  return (list(Xtr=Xtr, ytr=ytr, Xte=Xte, yte=yte, real_n=real_n))
}
