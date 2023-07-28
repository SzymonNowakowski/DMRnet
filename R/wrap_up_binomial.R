wrap_up_binomial <- function(mm, p, maxp, SS, fl, x.full, ord, n, levels.listed, mL, arguments) {
  maxl <- max(sapply(1:length(mm), function(i) length(mm[[i]]$loglik)))
  loglik <- sapply(1:length(mm), function(i) c(rep(-Inf, maxl - length(mm[[i]]$loglik)), mm[[i]]$loglik))

  out    <- prevent_merging_levels(clust_method = arguments$clust.method, result_matrix = loglik, as_vector = (maxl==1), mm = mm, SS = SS, factor_levels = fl, filler = -Inf)
  loglik <- out$result_matrix
  shift  <- out$shift

  out           <- calculate_model_groups(result_matrix = loglik, selection_function = which.max, p = p, maxp = maxp, shift = shift)
  idx           <- out$idx
  model_group   <- out$model_group
  model_index_within_group <- out$model_index_within_group

  be <- sapply(idx, function(i) {
    b_matrix<-mm[[model_group[i]]]$b;
    if (is.null(dim(b_matrix))) {
      b_matrix<-matrix(b_matrix);   #note this shouldn't be b_matrix<-t(as.matrix(b_matrix)) :   with other matrices that degenerated to HORIZONTAL vectors we want to have HORIZONTAL matrices
      #in this case, however, we want a VERTICAL matrix
      #in b_matrix[,model_index_within_group(i)] we take columns of b_matrix if b_matrix is a legitimate matrix
      #when it is degenerate (a vector), we want that taking the full first column (as model_index_within_group(i) == 1 in those cases) takes this whole vector
    }
    if (!is.na(model_index_within_group[i])) {
      return(part2beta_glm_help(b = b_matrix[,model_index_within_group[i]], S = SS[,model_group[i]], fl=fl))
    } else {
      return(rep(NA, 1+sum(fl-1)))
    }
  })

  legal_cols <- !is.na(apply(be, 2, sum))

  rownames(be) <- colnames(x.full)

  be <- be[ord,]  #reordering betas to reflect the original matrix X

  fit <- list(beta = be[,legal_cols], df = (length(idx):1)[legal_cols], loglik = loglik[cbind(idx[legal_cols], model_group[idx[legal_cols]])], n = n, levels.listed = levels.listed, lambda = mL$lambda, arguments=arguments, interc = TRUE)

  class(fit) = "DMR"
  return(fit)
}
