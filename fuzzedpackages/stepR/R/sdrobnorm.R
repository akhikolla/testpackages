sdrobnorm <- function(x, p = c(0.25, 0.75), lag = 1, supressWarningNA = FALSE, supressWarningResultNA = FALSE) {
  if (!is.numeric(x)) {
    stop("data vector 'x' must be a numeric vector")
  }
  
  if (!is.numeric(p) || length(p) != 2 || any(p < 0) || any(p > 1) || 
      isTRUE(all.equal(p[1], p[2], tolerance = 1e-6))) {
    stop("p must be a vector of two distinct probabilities, ",
         "i.e. a numeric vector with two distinct values between 0 and 1")
  }
  
  if (!is.numeric(lag) || length(lag) != 1) {
    stop("lag must be a single positive integer")
  }
  
  if (!is.integer(lag)) {
    lag <- as.integer(lag + 1e-6)
  }
  
  if (lag < 1) {
    stop("lag must be a single positive integer")
  }
  
  if (!is.logical(supressWarningNA) || length(supressWarningNA) != 1 || is.na(supressWarningNA)) {
    stop("supressWarningNA must be a single logical (not NA)")
  }
  
  if (!is.logical(supressWarningResultNA) || length(supressWarningResultNA) != 1 || 
      is.na(supressWarningResultNA)) {
    stop("supressWarningResultNA must be a single logical (not NA)")
  }
  
  if (any(is.na(x))) {
    x <- x[!is.na(x)]
    
    if (!supressWarningNA) {
      warning("the data vector 'x' contains NAs")
    }
  }
  
  if (length(x) < lag + 2L) {
    if (!supressWarningResultNA) {
      warning("result is NA, since the number of observations is too small")
    }
    return(NA_real_)
  }
  
  as.numeric(diff(quantile(diff(x, lag = lag), p)) / diff(qnorm(p)) / sqrt(2))
}
