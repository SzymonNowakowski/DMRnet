cv_indexation.mode_distribute <- function(X, y, nfolds, indexation.mode, model_function, ...) {

       X <- data.frame(X, check.names = TRUE, stringsAsFactors = TRUE)

       if (indexation.mode == "size") {
               out <- cv_size_indexed(X, y, nfolds, model_function, ...)
       } else {
               if (indexation.mode == "GIC") {
                       out <- cv_GIC_indexed(X, y, nfolds, model_function, ...)
               } else {
                       stop("Error: wrong indexation.mode, should be one of: GIC, size")
               }
       }


       class(out) <- "cv.DMR"

       return(out)
}
