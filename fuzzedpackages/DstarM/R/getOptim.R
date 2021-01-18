getOptim <- function(object) {
  if (!is.DstarM(object)) {
    stop("Input must be of class DstarM.", call. = FALSE)
  } else {
    name <- names(object)[1]
    if (name == "Bestvals") {
      return(object$GlobalOptimizer$optim$bestval)
    } else if (name == "r.hat") {
      return(sapply(object$GlobalOptimizer, function(x) x$optim$bestval))
    } else {
      warning("No optimizer information available for this object.")
      return(NULL)
    }
  }
}


