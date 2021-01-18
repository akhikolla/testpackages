BstrapTest <- function(y, bandwidth, nbandwidth = 30L, B = 500L, alpha = 0.05,
                          kernel = c("epanechnikov", "gaussian", "rectangular",
                                     "triangular", "biweight", "silverman")) {
  if (!is.numeric(y) || any(!is.finite(y))) {
    stop("observations 'y' must be a numeric vector containing only finite values")
  }
  
  n <- length(y)
  
  if (missing(bandwidth)) {
    if (!is.numeric(nbandwidth) || length(nbandwidth) != 1 || !is.finite(nbandwidth)) {
      stop("'nbandwidth' must be a single positive integer")
    }
    
    if (!is.integer(nbandwidth)) {
      nbandwidth <- as.integer(nbandwidth + 1e-6)
    }
    
    if (nbandwidth < 1L) {
      stop("'nbandwidth' must be a single positive integer")
    }
    
    bandwidth <- exp(seq(log(10 / n), log(0.25), length.out = nbandwidth))
  } else {
    if (!is.numeric(bandwidth) || any(!is.finite(bandwidth)) || any(bandwidth < 1 / n) || any(bandwidth > 0.5)) {
      stop("'bandwidth' must be a numeric vector containing only finite values between 1 / length(y) and 0.5")
    }
  }
  
  if (!is.function(kernel)) {
    kernel <- match.arg(kernel)
    
    kernel <- switch(kernel,
                     rectangular = function(x) 1 / 2,
                     triangular = function(x) 1 - abs(x),
                     epanechnikov = function(x) 3 / 4 * (1 - x^2),
                     biweight = function(x) 5 / 16 * (1 - x^2)^2,
                     gaussian = function(x) dnorm(x, 0, 1),
                     silverman = function(x) exp(-abs(x) / sqrt(2)) * sin(abs(x) / sqrt(2) + pi / 4) / 2,
                     stop("unknown kernel")
    )
  }
  
  if (!is.numeric(B) || length(B) != 1 || !is.finite(B)) {
    stop("'B' must be a single positive integer")
  }
  
  if (!is.integer(B)) {
    B <- as.integer(B + 1e-6)
  }
  
  if (B < 1L) {
    stop("'B' must be a single positive integer")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a probability, i.e., a single numeric between 0 and 1")
  }
  
  .BstrapTest(y = y, bandwidth = bandwidth, kernel = kernel, B = B, alpha = alpha)
}

.BstrapTest <- function(y, bandwidth, kernel, B, alpha, n = length(y)) {
  ret <- .estimateSingleCp(y = y, bandwidth = bandwidth, kernel = kernel, n = n)
  names(ret)[1] <- "piecewiseSignal"
  
  if (length(bandwidth) > 1) {
    cv <- numeric(length(bandwidth))
    for (i in seq_along(bandwidth)) {
      b <- as.integer(n * bandwidth[i] + 1e-12)
      cv[i] <- .CVtwosided(Y = y, K = kernel(+(1:b) / (n * bandwidth[i])))
    }
    ret$bandwidthSmooth <- bandwidth[which.min(cv)]
  } else {
    ret$bandwidthSmooth <- bandwidth
  }
  
  b <- as.integer(n * ret$bandwidthSmooth + 1e-12)
  ret$smoothSignal <- .kernelSmoothing(y, kernel((-b:b) / (n * ret$bandwidthSmooth)))
  
  etilde <- y - ret$piecewiseSignal
  ehat <- etilde - mean(etilde)
  
  Tstar <- numeric(B)
  for (b in 1:B) {
    eStar <- sample(ehat, length(y), replace = TRUE)
    
    yStar <- ret$smoothSignal + eStar
    
    Tstar[b] <- .estimateSingleCp(y = yStar, bandwidth = bandwidth, kernel = kernel, n = n)$size
  }
  
  ret$critVal <- quantile(abs(Tstar), 1 - alpha)
  ret$pValue <- mean(abs(ret$size) <= abs(Tstar))
  ret$outcome <- abs(ret$size) > as.numeric(ret$critVal)
  ret
}
  