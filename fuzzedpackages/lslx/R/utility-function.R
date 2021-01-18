pen <- function(start, penalty, set, weight) {
  if (missing(start)) {
    start <- NA
  }
  if (missing(penalty)) {
    penalty <- "default"
  } else {
    if (!(penalty %in% c("lasso", "ridge", "default"))) {
      stop("Argument 'penalty' in 'pen()' must be 'lasso', 'ridge', or 'default'.")      
    }
  }
  if (missing(set)) {
    set = 1
  } else {
    if (!is.numeric(set)) {
      stop("Argument 'set' in 'pen()' must be 1 or 2.")
    }
  }
  if (missing(weight)) {
    weight <- 1
  } else {
    if (!is.numeric(weight)) {
      stop("Argument 'weight' in 'pen()' must be numeric.")
    } else {
      if (weight < 0) {
        stop("Argument 'weight' in 'pen()' must be positive.")
      }
    }
  }
return(paste0("pen", "(", 
              "start=", start, ",", 
              "penalty=", penalty, ",",
              "set=", set, ",",
              "weight=", weight, ")"))
}