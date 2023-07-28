prevent_merging_levels <- function(clust_method, result_matrix, as_vector, mm, SS, factor_levels, filler) {

  if (as_vector)
    result_matrix <- t(as.matrix(result_matrix))      #making rss a horizontal one-row matrix

  shift <- rep(0, ncol(result_matrix))   #calculating shift based on number of Infs, because we are about to add Infs
  for (i in 1:ncol(result_matrix)) {
    shift[i] <- sum(is.infinite(result_matrix[, i]))
  }

  if (clust_method == "variable_selection")
    for (col in 1:ncol(result_matrix)) {
      b_matrix <- mm[[col]]$b
      if (is.null(dim(b_matrix))) {
        b_matrix<-matrix(b_matrix)
      }
      for (row in (shift[col]+1):nrow(result_matrix))
        if (is.finite(result_matrix[row, col])) {

          #analyse if this cell results from full factors only

          b_vector<-b_matrix[, row - shift[col]]

          S_vector <- SS[, col]
          pos_in_b <- 1

          for (i in seq_along(S_vector))
            if (S_vector[i]==1) {
              #check if positions related to this factor are all the same - either 0 or >0
              #and if >0 than all different to prevent merging levels
              b_fragment <- b_vector[(pos_in_b+1):(pos_in_b+factor_levels[i]-1)]
              b_zeros <- (b_fragment == 0)
              if (length(unique(b_zeros))!=1 | length(unique(b_fragment))!=(factor_levels[i]-1))  # mix of 0 and >0  OR not all different
                result_matrix[row, col] <- filler
              pos_in_b <- pos_in_b + factor_levels[i]-1
            }
        }
    }
  return(list(result_matrix = result_matrix, shift=shift))
}
