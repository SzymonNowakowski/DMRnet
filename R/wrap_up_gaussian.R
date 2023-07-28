wrap_up_gaussian <- function(mm, p, maxp, SS, fl, X, y, x.full, ord, n, levels.listed, mL, arguments) {
  maxl <- max(sapply(1:length(mm), function(i) length(mm[[i]]$rss)))
  rss <- sapply(1:length(mm), function(i) c(rep(Inf, maxl - length(mm[[i]]$rss)), mm[[i]]$rss))

  if (maxl == 1)
    rss <- t(as.matrix(rss))      #making rss a horizontal one-row matrix

  shift <- rep(0, ncol(rss))   #calculating shift based on number of Infs, because we are about to add Infs
  for (i in 1:ncol(rss)) {
    shift[i] <- sum(rss[, i]==Inf)
  }

  if (arguments$clust.method == "variable_selection")
    for (col in 1:ncol(rss)) {
      b_matrix <- mm[[col]]$b
      if (is.null(dim(b_matrix))) {
        b_matrix<-matrix(b_matrix)
      }
      for (row in (shift[col]+1):nrow(rss))
        if (is.finite(rss[row, col])) {

          #analyse if this cell results from full factors only

          b_vector<-b_matrix[, row - shift[col]]

          S_vector <- SS[, col]
          pos_in_b <- 1

          for (i in seq_along(S_vector))
            if (S_vector[i]==1) {
              #check if positions related to this factor are all the same - either 0 or >0
                                                      #and if >0 than all different to prevent merging levels
              b_fragment <- b_vector[(pos_in_b+1):(pos_in_b+fl[i]-1)]
              b_zeros <- (b_fragment == 0)
              if (length(unique(b_zeros))!=1 | length(unique(b_fragment))!=(fl-1))  # mix of 0 and >0  OR not all different
                rss[row, col] <- Inf
              pos_in_b <- pos_in_b + fl[i]-1
            }
        }
    }

  ind <- apply(rss, 1, which.min)  #in each row, which is the index of a model minimizing rss

  maxi <- min(p, maxp)
  if (length(ind) > maxi){
    idx <- (length(ind) - maxi):length(ind)   #but real model sizes are still length(idx):1
  } else {
    idx <- 1:length(ind)
  }

  #smallest models are last
  model_group <- ind    #we will simply call it differently (alias)
  model_index_within_group <- rep(0, length(ind))
  for (i in idx) {
    if (is.infinite(rss[i, model_group[i]])) {
      model_index_within_group[i] <- NA
    } else {
      model_index_within_group[i] <- i-shift[model_group[i]]
    }
  }

  be <- sapply(idx, function(i) {
    b_matrix<-mm[[model_group[i]]]$b;
    if (is.null(dim(b_matrix))) {
      b_matrix<-matrix(b_matrix);    #note this shouldn't be b_matrix<-t(as.matrix(b_matrix)) :   with other matrices that degenerated to HORIZONTAL vectors we want to have HORIZONTAL matrices
      #in this case, however, we want a VERTICAL matrix
      #in b_matrix[,model_index_within_group[i]] we take columns of b_matrix if b_matrix is a legitimate matrix
      #when it is degenerate (a vector), we want that taking the full first column (as model_index_within_group[i] == 1 in those cases) takes this whole vector
    }
    if (!is.na(model_index_within_group[i])) {
      return(part2beta_help(b = b_matrix[, model_index_within_group[i]], S = SS[, model_group[i]], X = X, y = y, fl=fl))
    } else {
      return(rep(NA, 1+sum(fl-1)))
    }
  })

  legal_cols <- !is.na(apply(be, 2, sum))

  rownames(be) <- colnames(x.full)

  be <- be[ord,]  #reordering betas to reflect the original matrix X

  fit <- list(beta = be[,legal_cols], df = (length(idx):1)[legal_cols], rss = rss[cbind(idx[legal_cols], ind[idx[legal_cols]])], n = n, levels.listed = levels.listed, lambda = mL$lambda, arguments = arguments, interc = TRUE)

  class(fit) = "DMR"
  return(fit)
}
