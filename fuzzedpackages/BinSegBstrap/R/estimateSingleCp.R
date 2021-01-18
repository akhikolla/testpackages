estimateSingleCp <- function(y, bandwidth, nbandwidth = 30L, 
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

    bandwidth <- exp(seq(log(2 / n), log(0.25), length.out = nbandwidth))
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
  
  .estimateSingleCp(y = y, bandwidth = bandwidth, kernel = kernel)
}

.estimateSingleCp <- function(y, bandwidth, kernel, n = length(y)) {
  if (length(bandwidth) > 1) {
    cv <- numeric(length(bandwidth))
    cp <- numeric(length(bandwidth))
    
    for (i in seq_along(bandwidth)) {
      b <- as.integer(n * bandwidth[i] + 1e-12)
      if (length(y) < 2L * b) {
        cv[i] <- Inf
        cp[i] <- NA
      } else {
        cp[i] <- as.integer(which.max(DmaxVec(y, b)) + b + 1e-12)
        cv[i] <- .CVonesided(Y = y[1:(cp[i] - 1)], K = kernel(+(1:b) / (n * bandwidth[i]))) + 
          .CVonesided(Y = rev(y[cp[i]:length(y)]), K = kernel(+(1:b) / (n * bandwidth[i])))
      }
    }

    if (all(is.na(cp))) {
      stop("all bandwidth values are too large")
    }
    
    indexBandwidth <- which.min(cv)
    bandwidth <- bandwidth[indexBandwidth]
    cp <- cp[indexBandwidth]
  } else {
    b <- as.integer(n * bandwidth + 1e-12)
    if (length(y) < 2L * b) {
      stop("bandwidth is too large")
    }
    cp <- as.integer(which.max(DmaxVec(y, b)) + b + 1e-12)
  }
  
  est <- numeric(length(y))
  b <- as.integer(n * bandwidth + 1e-12)
  est[1:(cp - 1)] <- .kernelSmoothing(y[1:(cp - 1)], kernel((-b:b) / (n * bandwidth)))
  est[cp:length(y)] <- .kernelSmoothing(y[cp:length(y)], kernel((-b:b) / (n * bandwidth)))
  size <- est[cp] - est[cp - 1]

  list(est = est, cp = cp, size = size, bandwidth = bandwidth)
}
  