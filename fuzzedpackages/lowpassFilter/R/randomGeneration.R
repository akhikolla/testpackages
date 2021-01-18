randomGeneration <- function(n, filter, signal = 0, noise = 1, oversampling = 100L, seed = n,
                             startTime = 0, truncated = TRUE) {
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n)) {
    stop("required number of observations 'n' must be a single positive integer")
  }
  
  if (!is.integer(n)) {
    n <- as.integer(n + 1e-6)
  }
  
  if (n < 1L) {
    stop("required number of observations 'n' must be a single positive integer")
  }

  if (!methods::is(filter, "lowpassFilter")) {
    stop("filter must be an object of class 'lowpassFilter'")
  }

  if (!is.numeric(startTime) || length(startTime) != 1 || !is.finite(startTime)) {
    stop("startTime must be a single finite numeric")
  }
  
  if (!is.logical(truncated) || length(truncated) != 1 || is.na(truncated)) {
    stop("truncated must be a single logical (not NA)")
  }

  if (!is.numeric(oversampling) || length(oversampling) != 1 || !is.finite(oversampling)) {
    stop("oversampling must be a single positive integer")
  }
  
  if (!is.integer(oversampling)) {
    oversampling <- as.integer(oversampling + 1e-6)
  }
  
  if (oversampling < 1L) {
    stop("oversampling must be a single positive integer")
  }

  if (methods::is(signal, "stepblock")) {
    time <- startTime + 1:n / filter$sr
    signal <- getConvolution(time, signal, filter, truncated)
  } else {
    if (!is.numeric(signal) || !all(is.finite(signal))) {
      stop("signal must be a finite numeric or an object of class 'stepblock'")
    }

    if (length(signal) != n && length(signal) != 1L) {
      stop("signal must be of length 1 or of length n if it is a numeric")
    }
  }

  if (methods::is(noise, "stepblock")) {
    time <- startTime + 
      (seq(1 - filter$len + 1 / oversampling, n, 1 / oversampling) - 1 / 2 / oversampling) / filter$sr
    
    ret <- as.numeric(rep(NA, length(time)))
    ret[time < noise$leftEnd[1]] <- noise$value[1]
    for (i in seq(along = noise$value)) {
      ret[time >= noise$leftEnd[i] & time < noise$rightEnd[i]] <- noise$value[i]
    }
    ret[time >= noise$rightEnd[i]] <- noise$value[i]
    noise <- ret
  } else {
    if (!is.numeric(noise) || !all(is.finite(noise))) {
      stop("noise must be a finite numeric or an object of class 'stepblock'")
    }
    
    if (length(noise) != (n + filter$len - 1L) * oversampling && length(noise) != 1L) {
      stop("noise must be of length 1 or of length ", (n + filter$len - 1L) * oversampling,
           " ((n + filter$len - 1L) * oversampling) if it is a numeric")
    }
  }
  
  if (seed == "no") {
    # set no seed
  } else {
    set.seed(seed)
  }

  diffFilter <- diff(filter$truncatedStepfun(seq(0, filter$len, 1 / oversampling) / filter$sr))
  eps <- stats::rnorm((n + filter$len - 1) * oversampling, 0, noise / sqrt(sum(diffFilter^2)))

  y <- signal + .convolveOversampling(eps, diffFilter, oversampling)
}

randomGenerationMA <- function(n, filter, signal = 0, noise = 1, seed = n, startTime = 0, truncated = TRUE) {
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n)) {
    stop("required number of observations 'n' must be a single positive integer")
  }
  
  if (!is.integer(n)) {
    n <- as.integer(n + 1e-6)
  }
  
  if (n < 1L) {
    stop("required number of observations 'n' must be a single positive integer")
  }
  
  if (!methods::is(filter, "lowpassFilter")) {
    stop("filter must be an object of class 'lowpassFilter'")
  }
  
  if (!is.numeric(startTime) || length(startTime) != 1 || !is.finite(startTime)) {
    stop("startTime must be a single finite numeric")
  }
  
  if (!is.logical(truncated) || length(truncated) != 1 || is.na(truncated)) {
    stop("truncated must be a single logical (not NA)")
  }
  
  if (methods::is(signal, "stepblock")) {
    time <- startTime + 1:n / filter$sr
    signal <- getConvolution(time, signal, filter, truncated)
  } else {
    if (!is.numeric(signal) || !all(is.finite(signal))) {
      stop("signal must be a finite numeric or an object of class 'stepblock'")
    }
    
    if (length(signal) != n && length(signal) != 1L) {
      stop("signal must be of length 1 or of length n if it is a numeric")
    }
  }
  
  if (!is.numeric(noise) || length(noise) != 1 || !is.finite(noise) || noise <= 0) {
    stop("noise must be a single postive finite numeric")
  }
  
  if (seed == "no") {
    # set no seed
  } else {
    set.seed(seed)
  }
  
  signal + .rand.genMDependent(n = n, filter = filter) * noise
}

.rand.genMDependent <- function(n, filter) {
  kern <- c(1, .computeMA(filter$acf))
  z <- stats::rnorm(n + length(kern) - 1, sd = 1)
  .convolve(z, kern) / sqrt(sum(kern^2))
}

.computeMA <- function(cov) {
  N <- length(cov) - 1
  alpha <- matrix(rep(0, N * N), nrow = N, ncol = N)
  vs <- rep(cov[1], N + 1)
  
  for(i in 1:N){
    alpha[i, 1:i] <- rep(0, i)
    alpha[i, i] <- cov[i + 1] / vs[1]
    if(i > 1) {
      for(k in 1:(i - 1)) {
        js <- 0:(k - 1)
        alpha[i, i - k] <- (cov[i - k + 1] - sum(alpha[i, i - js] * alpha[k, k - js] * vs[js + 1])) / vs[k + 1]
      }
    }
    js <- 0:(i - 1)
    vs[i + 1] <- vs[i + 1] - sum(alpha[i, i - js]^2 * vs[js + 1])
  }
  
  alpha[N, ]
}
