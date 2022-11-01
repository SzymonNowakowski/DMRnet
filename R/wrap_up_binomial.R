wrap_up_binomial <- function(mm, p, maxp, SS, fl, x.full, ord, n, levels.listed, mL, arguments) {
  maxl <- max(sapply(1:length(mm), function(i) length(mm[[i]]$loglik)))
  loglik <- sapply(1:length(mm), function(i) c(rep(-Inf, maxl - length(mm[[i]]$loglik)), mm[[i]]$loglik))

  if (maxl == 1)
    loglik <- t(as.matrix(loglik))      #making loglik a horizontal one-row matrix

  ind <- apply(loglik, 1, which.max)  #in each row, which is the index of a model maximizing loglik

  maxi <- min(p, maxp)
  if (length(ind) > maxi){
    idx <- (length(ind) - maxi):length(ind)
  } else{
    idx <- 1:length(ind)
  }

  shift<- function(i) {sum(loglik[, i]==-Inf)}
  model_group <- function(i) {ind[i]}
  model_index_within_group<- function(i) {i-shift(model_group(i))}

  be <- sapply(idx, function(i) {
    b_matrix<-mm[[model_group(i)]]$b;
    if (is.null(dim(b_matrix))) {
      b_matrix<-matrix(b_matrix);   #note this shouldn't be b_matrix<-t(as.matrix(b_matrix)) :   with other matrices that degenerated to HORIZONTAL vectors we want to have HORIZONTAL matrices
      #in this case, however, we want a VERTICAL matrix
      #in b_matrix[,model_index_within_group(i)] we take columns of b_matrix if b_matrix is a legitimate matrix
      #when it is degenerate (a vector), we want that taking the full first column (as model_index_within_group(i) == 1 in those cases) takes this whole vector
    }
    part2beta_glm_help(b = b_matrix[,model_index_within_group(i)], S = SS[,model_group(i)], fl=fl)
  })

  rownames(be) <- colnames(x.full)

  be <- be[ord,]  #reordering betas to reflect the original matrix X

  fit <- list(beta = be, df = length(idx):1, loglik = loglik[cbind(idx, ind[idx])], n = n, levels.listed = levels.listed, lambda = mL$lambda, arguments=arguments, interc = TRUE)
  class(fit) = "DMR"
  return(fit)
}
