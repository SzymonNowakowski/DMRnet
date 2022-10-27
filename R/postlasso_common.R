postlasso_common <- function(lambdas, n, betas) {

  #called by DMRnet, glamer and SOSnet

  dfy <- apply(betas, 2, function(x) sum(x!=0))
  kt <- 1:length(lambdas)
  if (dfy[1] != 1L) { ##checking that maybe the largest lambda is NOT the-only-intercept lambda
    # cases like this have been observed in high dimensional simulations (https://github.com/SzymonNowakowski/DMRnet/issues/39)
    #if it is so, we must add an artificial the-only-intercept lambda
    # it is a strictly technical add-on, i.e. it results in SS being calculated later that it includes a full model as the first component
    #so even if interc=FALSE in SOSnet we proceed like this, the first row will be removed from bb anyway
                  # just outside of this function in SOSnet


    betas <- cbind(c(1, rep(0, nrow(betas)-1)), betas)
    colnames(betas) <- c("artificial_only_intercept_lambda", colnames(betas)[-1])

    dfy <- apply(betas, 2, function(x) sum(x!=0))   #recalculate the dfy & kt
    kt <- 1:(length(lambdas)+1)
  }

  lambdas_with_nonzero_beta_number_too_large <- which(dfy >= n)  #(1) removing predictor sets with more predictors than matrix rows
  if (length(lambdas_with_nonzero_beta_number_too_large) > 0){
    kt <- kt[-lambdas_with_nonzero_beta_number_too_large]  #removing them from lambda indices
    dfy <- dfy[-lambdas_with_nonzero_beta_number_too_large]
  }
  lambdas_with_no_betas <- which(dfy == 0)     #(2) removing predictor sets with 0 predictors
  if(length(lambdas_with_no_betas) > 0){
    kt <- kt[-lambdas_with_no_betas]  #removing them from lambda indices
  }
  bb <- as.matrix(abs(betas[, kt]))  #bb is a matrix listing beta values (rows) respective to the net of lambda values (cols)
  bb_predictor_sets <- ifelse(bb > 0, 1, 0)          #bb_predictor_sets is a matrix listing predictor sets (0 or 1 for each predictor)  (rows) respective to the net of lambda values (cols)
  ii <- duplicated(t(bb_predictor_sets))    #detecting duplicated predictor sets


  bb <- bb[, ii == FALSE, drop = FALSE]    #(3) removing duplicated predictor sets

  return(bb)
}
