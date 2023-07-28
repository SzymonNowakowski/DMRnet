wrap_up_gaussian <- function(mm, p, maxp, SS, fl, X, y, x.full, ord, n, levels.listed, mL, arguments) {
  maxl <- max(sapply(1:length(mm), function(i) length(mm[[i]]$rss)))
  rss <- sapply(1:length(mm), function(i) c(rep(Inf, maxl - length(mm[[i]]$rss)), mm[[i]]$rss))

  out   <- prevent_merging_levels(clust_method = arguments$clust.method, result_matrix = rss, as_vector = (maxl==1), mm = mm, SS = SS, factor_levels = fl, filler = Inf)
  rss   <- out$result_matrix
  shift <- out$shift

  out           <- calculate_model_groups(result_matrix = rss, selection_function = which.min, p = p, maxp = maxp, shift = shift)
  idx           <- out$idx
  model_group   <- out$model_group
  model_index_within_group <- out$model_index_within_group

  be <- sapply(idx, function(i) {
    b_matrix<-as.matrix(mm[[model_group[i]]]$b)
      #note this shouldn't be b_matrix<-t(as.matrix(b_matrix)) :   with other matrices that degenerated to HORIZONTAL vectors we want to have HORIZONTAL matrices
      #in this case, however, we want a VERTICAL matrix
      #in b_matrix[,model_index_within_group[i]] we take columns of b_matrix if b_matrix is a legitimate matrix
      #when it is degenerate (a vector), we want that taking the full first column (as model_index_within_group[i] == 1 in those cases) takes this whole vector

    if (!is.na(model_index_within_group[i])) {
      return(part2beta_help(b = b_matrix[, model_index_within_group[i]], S = SS[, model_group[i]], X = X, y = y, fl=fl))
    } else {
      return(rep(NA, 1+sum(fl-1)))
    }
  })

  legal_cols <- !is.na(apply(be, 2, sum))

  rownames(be) <- colnames(x.full)

  be <- as.matrix(be[ord,])  #reordering betas to reflect the original matrix X, making it a matrix just in case it is one column

  fit <- list(beta = as.matrix(be[,legal_cols]), df = (length(idx):1)[legal_cols], rss = rss[cbind(idx[legal_cols], model_group[idx[legal_cols]])], n = n, levels.listed = levels.listed, lambda = mL$lambda, arguments = arguments, interc = TRUE)

  class(fit) = "DMR"
  return(fit)
}
