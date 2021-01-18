BinSegBstrap <- function(y, bandwidth, nbandwidth = 30L, B = 500, alpha = 0.05,
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

  ret <- .BinSegBstrap(y = y, kernel = kernel, bandwidth = bandwidth, B = B, alpha = alpha)
  
  cps <- c(1, ret, n + 1)
  
  if (length(bandwidth) > 1) {
    cv <- numeric(length(bandwidth))
    for (i in seq_along(bandwidth)) {
      b <- as.integer(n * bandwidth[i] + 1e-12)
      
      for (k in 1:(length(cps) - 1)) {
        cv[i] <- cv[i] + .CVtwosided(Y = y[cps[k]:(cps[k + 1] - 1)], K = kernel(+(1:b) / (n * bandwidth[i])))
      }
    }
    bandwidth <- bandwidth[which.min(cv)]
  }
  
  est <- numeric(n)
  b <- as.integer(n * bandwidth + 1e-12)

  for (i in 1:(length(cps) - 1)) {
    est[cps[i]:(cps[i + 1] - 1)] <- .kernelSmoothing(y[cps[i]:(cps[i + 1] - 1)], kernel((-b:b) / (n * bandwidth)))
  }
  
  list(est = est, cps = ret, bandwidth = bandwidth)
}

.BinSegBstrap <- function(y, kernel, bandwidth, B, alpha, n = length(y), s = 1, e = length(y)) {
  if (e - s <= 2 * as.integer(min(bandwidth) * n + 1e-12)) {
    return(integer(0))
  } else {
    ret <- .BstrapTest(y = y[s:e], kernel = kernel, bandwidth = bandwidth, B = B, alpha = alpha, n = n) 
    
    if (ret$outcome) {
      loc <- ret$cp + s - 1L
      
      cpLeft <- .BinSegBstrap(y = y, kernel = kernel, bandwidth = bandwidth, B = B, alpha = alpha,
                              n = n, s = s, e = loc - 1)
      cpRight <- .BinSegBstrap(y = y, kernel = kernel, bandwidth = bandwidth, B = B, alpha = alpha,
                               n = n, s = loc, e = e)
    } else {
      return(integer(0))
    }
  }
  c(cpLeft, loc, cpRight)
}
