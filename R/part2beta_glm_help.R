part2beta_glamer_glm_help <- function(b, S, fl, valid = TRUE) {
  if (!valid) {
    return(rep(0, sum(fl-1)+1))
  }

  if (sum(S == 0) > 0){
    b1 <- b[1]
    b <- b[-1]
    for (i in 1:length(S)){
      if(S[i] == 1) {
        b1 <- c(b1, b[1:(fl[i] - 1)])
        b <- b[-c(1:(fl[i] - 1))]
      } else{
        b1 <- c(b1, rep(0, (fl[i] - 1)))
      }

    }
  } else{
    b1 <- b
  }
  return(b1)
}
