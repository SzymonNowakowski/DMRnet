wrap_up_gaussian <- function(mm, p, maxp, SS, fl, X, y, x.full, ord, n, levels.listed, mL, arguments) {
  maxl <- max(sapply(1:length(mm), function(i) length(mm[[i]]$rss)))
  rss <- sapply(1:length(mm), function(i) c(rep(Inf, maxl - length(mm[[i]]$rss)), mm[[i]]$rss))

  if (maxl == 1)
    rss <- t(as.matrix(rss))      #making rss a horizontal one-row matrix

  ind <- apply(rss, 1, which.min)  #in each row, which is the index of a model minimizing rss

  maxi <- min(p, maxp)
  if (length(ind) > maxi){
    idx <- (length(ind) - maxi):length(ind)   #but real model sizes are still length(idx):1
  } else {
    idx <- 1:length(ind)
  }

  #smallest models are last
  shift<- function(i) {sum(rss[, i]==Inf)}
  model_group <- function(i) {ind[i]}
  model_index_within_group<- function(i) {i-shift(model_group(i))}

  be <- sapply(idx, function(i) {
    b_matrix<-mm[[model_group(i)]]$b;
    if (is.null(dim(b_matrix))) {
      b_matrix<-matrix(b_matrix);    #note this shouldn't be b_matrix<-t(as.matrix(b_matrix)) :   with other matrices that degenerated to HORIZONTAL vectors we want to have HORIZONTAL matrices
      #in this case, however, we want a VERTICAL matrix
      #in b_matrix[,model_index_within_group(i)] we take columns of b_matrix if b_matrix is a legitimate matrix
      #when it is degenerate (a vector), we want that taking the full first column (as model_index_within_group(i) == 1 in those cases) takes this whole vector
    }
    part2beta_help(b = b_matrix[, model_index_within_group(i)], S = SS[, model_group(i)], X = X, y = y, fl=fl)
  })

  rownames(be) <- colnames(x.full)
  if(length(ord) > 0){
    ind1 <- c(1:p)
    ind1[-ord] = 1:(p - length(ord))
    ind1[ord] = (p - length(ord) + 1):p
    be = be[ind1,]
  }

  fit <- list(beta = be, df = length(idx):1, rss = rss[cbind(idx, ind[idx])], n = n, levels.listed = levels.listed, lambda = mL$lambda, arguments = arguments, interc = TRUE)

  class(fit) = "DMR"
  return(fit)
}
