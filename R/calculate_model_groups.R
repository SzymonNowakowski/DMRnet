calculate_model_groups <- function(result_matrix, selection_function, p, maxp, shift) {
  ind <- apply(result_matrix, 1, selection_function)  #in each row, which is the index of a model maximizing loglik
                                                                  #           or index minimizing RSS

  maxi <- min(p, maxp)
  if (length(ind) > maxi){
    idx <- (length(ind) - maxi):length(ind)
  } else {
    idx <- 1:length(ind)
  }

  #smallest models are last

  model_index_within_group <- rep(0, length(ind))
  for (i in idx) {
    if (is.infinite(result_matrix[i, ind[i]])) {
      model_index_within_group[i] <- NA
    } else {
      model_index_within_group[i] <- i-shift[ind[i]]
    }
  }
  return(list(idx = idx, model_group = ind, model_index_within_group = model_index_within_group))
}
