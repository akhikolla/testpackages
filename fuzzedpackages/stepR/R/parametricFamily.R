.parametricFamily <- function(family, y = NULL, n, nq = 2L^(as.integer(log2(n) + 1e-12) + 1L) - 1L, ...) {
  if (is.null(family)) {
    family <- "gauss"
  }
  family <- match.arg(family, c("gauss", "mDependentPS", "jsmurf", "jsmurfPS", "jsmurfLR",
                                "hsmuce", "hjsmurf", "hjsmurfSPS", "hjsmurfLR",
                                "2Param", "LR"))
  
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n)) {
    stop("number of observations 'n' must be a single positive integer")
  }
  
  if (!is.integer(n)) {
    n <- as.integer(n + 1e-6)
  }
  
  if (n < 1L) {
    stop("number of observations 'n' must be a single positive integer")
  }
  
  if (!is.numeric(nq) || length(nq) != 1 || !is.finite(nq)) {
    stop("nq must be a single integer greather or equal than the number of observations 'n'")
  }
  
  if (!is.integer(nq)) {
    nq <- as.integer(nq + 1e-6)
  }
  
  if (nq < n) {
    stop("nq must be a single finite integer greather or equal than the number of observations 'n'")
  }
  
  if (is.null(y)) {
    data <- list(family = family, n = n, nq = nq)
  } else {
    data <- list(y = NULL, family = family, n = n, nq = nq)
  }
  
  switch(family,
         "gauss" = .familyGauss(data = data, y = y, ...),
         "mDependentPS" = .familyMDependentPS(data = data, y = y, ...),
         "jsmurf" = .familyJsmurf(data = data, y = y, ...),
         "jsmurfPS" = .familyJsmurfPS(data = data, y = y, ...),
         "jsmurfLR" = .familyJsmurfLR(data = data, y = y, ...),
         "hsmuce" = .familyHsmuce(data = data, y = y, ...),
         "hjsmurf" = .familyHjsmurf(data = data, y = y, ...),
         "hjsmurfSPS" = .familyHjsmurfSPS(data = data, y = y, ...),
         "hjsmurfLR" = .familyHjsmurfLR(data = data, y = y, ...),
         "2Param" = .family2Param(data = data, y = y, ...),
         "LR" = .familyLR(data = data, y = y, ...),
         stop("unknown family"))
}

.familyGauss <- function(data, y, sd = NULL) {
  if (!is.null(y)) {
    if (!is.double(y) || any(!is.finite(y))) {
      stop("for family 'gauss' y must be a double vector containing only finite values")
    }
    
    if (length(y) != data$n) {
      stop("y must be of length n")
    }
    
    data$y <- y
  }
  
  if (is.null(sd)) {
    if (is.null(y)) {
      data$sd <- 1
    } else {
      data$sd <- sdrobnorm(y, supressWarningResultNA = TRUE)
      if (is.na(data$sd)) {
        stop("number of observations is too small for computing 'sd' by 'sdrobnorm'")
      }
    }
  } else {
    if (!is.numeric(sd) || length(sd) != 1 || !is.finite(sd) || sd <= 0) {
      stop("sd must be a single positive finite numeric")
    }
    
    data$sd <- as.numeric(sd)
  }
  data$type <- 0L
  data$argumentsList <- list(sd = data$sd)
  data$possibleLengths <- 1:data$n
  data$defaultLengths <- 1:data$n
  if (data$n <= 1e3L) {
    data$defaultIntervalSystem <- "all"
  } else {
    data$defaultIntervalSystem <- "dyaLen"
  }
  data$defaultPenalty <- "sqrt"
  
  data$key <- digest::sha1(list("gauss"), digits = 6)
  data$penaltyShift <- 0
  data$rand.gen <- function(data) rnorm(data$n)
  data$addData <- function(rand.gen, data) {
    data$y <- rand.gen(data)
    
    if (!is.numeric(data$y) || length(data$y) != data$n || any(!is.finite(data$y))) {
      stop("for parametric family 'gauss' rand.gen must return a finite numeric vector ",
           "of length equal to the number of observations")
    }
    data
  }
  data$generateSignal <- function(data, intervalSystem) {
    0
  }
  data$MC <- function(data, n) {
    data$n <- n
    data$nq <- n
    data$possibleLengths <- 1:n
    data$defaultLengths <- 1:n
    if (n <= 1e3L) {
      data$defaultIntervalSystem <- "all"
    } else {
      data$defaultIntervalSystem <- "dyaLen"
    }
    data$sd <- 1
    data$argumentsList <- list(sd = data$sd)
    data
  }
  data$save <- TRUE

  data
}

.familyMDependentPS <- function(data, y, covariances = NULL, correlations = NULL, filter = NULL, sd = NULL) {
  if (!is.null(y)) {
    if (!is.double(y) || !all(is.finite(y))) {
      stop("for family 'mDependentPS' y must be a double vector containing only finite values")
    }
    
    if (length(y) != data$n) {
      stop("y must be of length n")
    }
    
    data$y <- y
  }
  
  if (is.null(covariances)) {
    if (is.null(correlations)) {
      if (is.null(filter)) {
        stop("covariances, correlations or filter must be given for family 'mDependentPS'")
      } else {
        if (!is(filter, "lowpassFilter")) {
          stop("filter must be an object of class 'lowpassFilter'")
        }
        
        correlations <- filter$acf
      }
    } else {
      if (!is.numeric(correlations) || !all(is.finite(correlations))) {
        stop("correlations must be a finite numeric vector")
      }
      
      if (correlations[1] != 1 || any(abs(correlations[-1]) > 1) ||
          correlations[length(correlations)] == 0) {
        stop("correlations must be a correlation vector, ",
             "i.e. the first element must be 1, the absolute value of every other element must be ",
             "smaller or equal and the last element should not be zero")
      }
    }
    
    if (is.null(sd)) {
      if (is.null(y)) {
        sd <- 1
      } else {
        sd <- sdrobnorm(y, lag = length(correlations), supressWarningResultNA = TRUE)
        if (is.na(sd)) {
          stop("number of observations is too small for computing 'sd' by 'sdrobnorm'")
        }
      }
    } else {
      if (!is.numeric(sd) || length(sd) != 1 || !is.finite(sd) || sd <= 0) {
        stop("sd must be a single positive finite numeric")
      }
      
      sd <- as.numeric(sd)
    }
    
    covariances <- sd^2 * correlations
  } else {
    if (!is.numeric(covariances)  || !all(is.finite(covariances))) {
      stop("covariances has to be a finite numeric vector")
    }
    
    if (covariances[1] <= 0 || any(abs(covariances[-1]) > covariances[1]) ||
        covariances[length(covariances)] == 0) {
      stop("covariances has to be a covariance vector, ",
           "i.e. the first element has to be positive, the absolute value of every other element has to be ",
           "smaller or equal and the last element should not be zero")
    }
  }
  data$type <- 10L
  data$covariances <- covariances
  data$argumentsList <- list(covariances = data$covariances)
  data$possibleLengths <- 1:data$n
  data$defaultLengths <- 1:data$n
  data$defaultIntervalSystem <- "dyaLen"
  data$defaultPenalty <- "sqrt"
  data$ma <- .computeMA(data$covariances)
  data$penaltyShift <- 0
  data$rand.gen <- .rand.genMDependent
  data$addData <- function(rand.gen, data) {
    data$y <- rand.gen(data)
    
    if (!is.numeric(data$y) || length(data$y) != data$n || any(!is.finite(data$y))) {
      stop("for parametric family 'mDependentPS' rand.gen must return a finite numeric vector ",
           "of length equal to the number of observations")
    }
    data
  }
  data$key <- digest::sha1(list("mDependentPS", data$ma), digits = 6)
  data$generateSignal <- function(data, intervalSystem) {
    0
  }
  data$MC <- function(data, n) {
    data$n <- n
    data$nq <- n
    data$possibleLengths <- 1:n
    data$defaultLengths <- 1:n
    data$covariances <- data$covariances / data$covariances[1]
    data$argumentsList <- list(covariances = data$covariances)
    data
  }
  data$save <- TRUE

  data
}

.familyJsmurf <- function(data, y, sd = NULL, filter = NULL) {
  if (!is.null(y)) {
    if (!is.double(y) || any(!is.finite(y))) {
      stop("for family 'jsmurf' y must be a double vector containing only finite values")
    }
    
    if (length(y) != data$n) {
      stop("y must be of length n")
    }
    
    data$y <- y
  }
  
  if (is.null(filter)) {
    stop("filter must be given for family 'jsmurf'")
  } 
  
  if (!is(filter, "lowpassFilter")) {
    stop("filter must be an object of class 'lowpassFilter'")
  }
  
  if (is.null(sd)) {
    if (is.null(y)) {
      data$sd <- 1
    } else {
      data$sd <- sdrobnorm(y, lag = filter$len + 1L, supressWarningResultNA = TRUE)
      if (is.na(data$sd)) {
        stop("number of observations is too small to compute 'sd' by 'sdrobnorm'")
      }
    }
  } else {
    if (!is.numeric(sd) || length(sd) != 1 || !is.finite(sd) || sd <= 0) {
      stop("sd must be a single positive finite numeric")
    }
    
    data$sd <- as.numeric(sd)
  }
  
  data$type <- 11L
  data$filter <- filter
  data$argumentsList <- list(sd = data$sd, filterLength = data$filter$len)
  data$possibleLengths <- (data$filter$len + 1L):data$n
  data$defaultLengths <- (data$filter$len + 1L):data$n
  data$defaultIntervalSystem <- "dyaLen"
  data$defaultPenalty <- "sqrt"
  data$penaltyShift <- -data$filter$len
  
  data$ma <- .computeMA(data$filter$acf)
  data$key <- digest::sha1(list("jsmurf", data$ma), digits = 6)
  data$rand.gen <- .rand.genMDependent
  data$addData <- function(rand.gen, data) {
    data$y <- rand.gen(data)
    
    if (!is.numeric(data$y) || length(data$y) != data$n || any(!is.finite(data$y))) {
      stop("for parametric family 'jsmurf' rand.gen must return a finite numeric vector ",
           "of length equal to the number of observations")
    }
    data
  }
  data$generateSignal <- function(data, intervalSystem) {
    0
  }
  data$MC <- function(data, n) {
    data$n <- n
    data$nq <- n
    data$possibleLengths <- (data$filter$len + 1L):data$n
    data$defaultLengths <- (data$filter$len + 1L):data$n
    data$sd <- 1
    data$argumentsList <- list(sd = data$sd, filterLength = filter$len)
    data
  }
  data$save <- TRUE
  
  data
}

.familyJsmurfPS <- function(data, y, sd = NULL, filter = NULL) {
  if (!is.null(y)) {
    if (!is.double(y) || any(!is.finite(y))) {
      stop("for family 'jsmurf' y must be a double vector containing only finite values")
    }
    
    if (length(y) != data$n) {
      stop("y must be of length n")
    }
    
    data$y <- y
  }
  
  if (is.null(filter)) {
    stop("filter must be given for family 'jsmurfPS'")
  } 
  
  if (!is(filter, "lowpassFilter")) {
    stop("filter must be an object of class 'lowpassFilter'")
  }
  
  if (is.null(sd)) {
    if (is.null(y)) {
      data$sd <- 1
    } else {
      data$sd <- sdrobnorm(y, lag = filter$len + 1L, supressWarningResultNA = TRUE)
      if (is.na(data$sd)) {
        stop("number of observations is too small to compute 'sd' by 'sdrobnorm'")
      }
    }
  } else {
    if (!is.numeric(sd) || length(sd) != 1 || !is.finite(sd) || sd <= 0) {
      stop("sd must be a single positive finite numeric")
    }
    
    data$sd <- as.numeric(sd)
  }
  
  data$type <- 12L
  data$filter <- filter
  data$covariances <- data$sd * data$filter$acf
  data$argumentsList <- list(covariances = data$covariances, filterLength = data$filter$len)
  data$possibleLengths <- (data$filter$len + 1L):data$n
  data$defaultLengths <- (data$filter$len + 1L):data$n
  data$defaultIntervalSystem <- "dyaLen"
  data$defaultPenalty <- "sqrt"
  data$penaltyShift <- -data$filter$len
  data$ma <- .computeMA(data$filter$acf)
  data$key <- digest::sha1(list("jsmurfPS", data$ma), digits = 6)
  data$rand.gen <- .rand.genMDependent
  data$addData <- function(rand.gen, data) {
    data$y <- rand.gen(data)
    
    if (!is.numeric(data$y) || length(data$y) != data$n || any(!is.finite(data$y))) {
      stop("for parametric family 'jsmurfPS' rand.gen must return a finite numeric vector ",
           "of length equal to the number of observations")
    }
    data
  }
  data$generateSignal <- function(data, intervalSystem) {
    0
  }
  data$MC <- function(data, n) {
    data$n <- n
    data$nq <- n
    data$possibleLengths <- (data$filter$len + 1L):data$n
    data$defaultLengths <- (data$filter$len + 1L):data$n
    data$covariances <- data$covariances / data$covariances[1]
    data$argumentsList <- list(covariances = data$covariances, filterLength = filter$len)
    data
  }
  data$save <- TRUE
  
  data
}

.familyJsmurfLR <- function(data, y, sd = NULL, filter = NULL) {
  if (!is.null(y)) {
    if (!is.double(y) || any(!is.finite(y))) {
      stop("for family 'jsmurf' y must be a double vector containing only finite values")
    }
    
    if (length(y) != data$n) {
      stop("y must be of length n")
    }
    
    data$y <- y
  }
  
  if (is.null(filter)) {
    stop("filter must be given for family 'jsmurfLR'")
  } 
  
  if (!is(filter, "lowpassFilter")) {
    stop("filter must be an object of class 'lowpassFilter'")
  }
  
  if (is.null(sd)) {
    if (is.null(y)) {
      data$sd <- 1
    } else {
      data$sd <- sdrobnorm(y, lag = filter$len + 1L, supressWarningResultNA = TRUE)
      if (is.na(data$sd)) {
        stop("number of observations is too small to compute 'sd' by 'sdrobnorm'")
      }
    }
  } else {
    if (!is.numeric(sd) || length(sd) != 1 || !is.finite(sd) || sd <= 0) {
      stop("sd must be a single positive finite numeric")
    }
    
    data$sd <- as.numeric(sd)
  }
  
  data$type <- 13L
  data$filter <- filter
  data$covariances <- data$sd * data$filter$acf
  data$argumentsList <- list(covariances = data$covariances, filterLength = data$filter$len)
  data$possibleLengths <- (data$filter$len + 1L):data$n
  data$defaultLengths <- (data$filter$len + 1L):data$n
  data$defaultIntervalSystem <- "dyaLen"
  data$defaultPenalty <- "sqrt"
  data$penaltyShift <- -data$filter$len
  data$ma <- .computeMA(data$filter$acf)
  data$key <- digest::sha1(list("jsmurfLR", data$ma), digits = 6)
  data$rand.gen <- .rand.genMDependent
  data$addData <- function(rand.gen, data) {
    data$y <- rand.gen(data)
    
    if (!is.numeric(data$y) || length(data$y) != data$n || any(!is.finite(data$y))) {
      stop("for parametric family 'jsmurfLR' rand.gen must return a finite numeric vector ",
           "of length equal to the number of observations")
    }
    data
  }
  data$generateSignal <- function(data, intervalSystem) {
    0
  }
  data$MC <- function(data, n) {
    data$n <- n
    data$nq <- n
    data$possibleLengths <- (data$filter$len + 1L):data$n
    data$defaultLengths <- (data$filter$len + 1L):data$n
    data$covariances <- data$covariances / data$covariances[1]
    data$argumentsList <- list(covariances = data$covariances, filterLength = filter$len)
    data
  }
  data$save <- TRUE
  
  data
}

.familyHsmuce <- function(data, y) {
  if (!is.null(y)) {
    if (!is.double(y) || !all(is.finite(y))) {
      stop("for family 'hsmuce' y must be a double vector containing only finite values")
    }
    
    if (length(y) != data$n) {
      stop("y must be of length n")
    }
    
    data$y <- y
  }
  
  data$type <- 20L
  data$argumentsList <- NULL
  data$possibleLengths <- 2:data$n
  data$defaultLengths <- 2:data$n
  data$defaultIntervalSystem <- "dyaPar"
  data$defaultPenalty <- "weights"
  data$rand.gen <- function(data) rnorm(data$n)
  data$penaltyShift <- 0
  data$addData <- function(rand.gen, data) {
    data$y <- rand.gen(data)
    
    if (!is.numeric(data$y) || length(data$y) != data$n || any(!is.finite(data$y))) {
      stop("for parametric family 'hsmuce' rand.gen must return a finite numeric vector ",
           "of length equal to the number of observations")
    }
    data
  }
  data$key <- digest::sha1(list("hsmuce"), digits = 6)
  data$generateSignal <- function(data, intervalSystem) {
    0
  }
  data$MC <- function(data, n) {
    data$n <- n
    data$nq <- n
    data$possibleLengths <- 2:n
    data$defaultLengths <- 2:n
    data
  }
  data$save <- TRUE

  data
}

.familyHjsmurf <- function(data, y, filter = NULL) {
  if (!is.null(y)) {
    if (!is.double(y) || any(!is.finite(y))) {
      stop("for family 'jsmurf' y must be a double vector containing only finite values")
    }
    
    if (length(y) != data$n) {
      stop("y must be of length n")
    }
    
    data$y <- y
  }
  
  if (is.null(filter)) {
    stop("filter must be given for family 'hjsmurf'")
  } 
  
  if (!is(filter, "lowpassFilter")) {
    stop("filter must be an object of class 'lowpassFilter'")
  }
  
  data$type <- 21L
  data$filter <- filter
  data$argumentsList <- list(filterLength = filter$len)
  data$possibleLengths <- (data$filter$len + 2L):data$n
  data$defaultLengths <- (data$filter$len + 2L):data$n
  data$defaultIntervalSystem <- "dyaLen"
  data$defaultPenalty <- "weights"
  data$penaltyShift <- -data$filter$len
  
  data$ma <- .computeMA(data$filter$acf)
  data$key <- digest::sha1(list("hjsmurf", data$ma), digits = 6)
  data$rand.gen <- .rand.genMDependent
  data$addData <- function(rand.gen, data) {
    data$y <- rand.gen(data)
    
    if (!is.numeric(data$y) || length(data$y) != data$n || any(!is.finite(data$y))) {
      stop("for parametric family 'hjsmurf' rand.gen must return a finite numeric vector ",
           "of length equal to the number of observations")
    }
    data
  }
  data$generateSignal <- function(data, intervalSystem) {
    0
  }
  data$MC <- function(data, n) {
    data$n <- n
    data$nq <- n
    data$possibleLengths <- (data$filter$len + 2L):data$n
    data$defaultLengths <- (data$filter$len + 2L):data$n
    data
  }
  data$save <- TRUE
  
  data
}

.familyHjsmurfSPS <- function(data, y, filter = NULL) {
  if (!is.null(y)) {
    if (!is.double(y) || any(!is.finite(y))) {
      stop("for family 'jsmurf' y must be a double vector containing only finite values")
    }
    
    if (length(y) != data$n) {
      stop("y must be of length n")
    }
    
    data$y <- y
  }
  
  if (is.null(filter)) {
    stop("filter must be given for family 'hjsmurfSPS'")
  } 
  
  if (!is(filter, "lowpassFilter")) {
    stop("filter must be an object of class 'lowpassFilter'")
  }
  
  data$type <- 22L
  data$filter <- filter
  data$correlations <- data$filter$acf
  data$argumentsList <- list(correlations = data$correlations, filterLength = filter$len)
  data$possibleLengths <- (data$filter$len + 2L):data$n
  data$defaultLengths <- (data$filter$len + 2L):data$n
  data$defaultIntervalSystem <- "dyaLen"
  data$defaultPenalty <- "weights"
  data$penaltyShift <- -data$filter$len
  
  data$ma <- .computeMA(data$filter$acf)
  data$key <- digest::sha1(list("hjsmurfSPS", data$ma), digits = 6)
  data$rand.gen <- .rand.genMDependent
  data$addData <- function(rand.gen, data) {
    data$y <- rand.gen(data)
    
    if (!is.numeric(data$y) || length(data$y) != data$n || any(!is.finite(data$y))) {
      stop("for parametric family 'hjsmurfPS' rand.gen must return a finite numeric vector ",
           "of length equal to the number of observations")
    }
    data
  }
  data$generateSignal <- function(data, intervalSystem) {
    0
  }
  data$MC <- function(data, n) {
    data$n <- n
    data$nq <- n
    data$possibleLengths <- (data$filter$len + 2L):data$n
    data$defaultLengths <- (data$filter$len + 2L):data$n
    data
  }
  data$save <- TRUE
  
  data
}

.familyHjsmurfLR <- function(data, y, filter = NULL) {
  if (!is.null(y)) {
    if (!is.double(y) || any(!is.finite(y))) {
      stop("for family 'jsmurf' y must be a double vector containing only finite values")
    }
    
    if (length(y) != data$n) {
      stop("y must be of length n")
    }
    
    data$y <- y
  }
  
  if (is.null(filter)) {
    stop("filter must be given for family 'hjsmurfLR'")
  } 
  
  if (!is(filter, "lowpassFilter")) {
    stop("filter must be an object of class 'lowpassFilter'")
  }
  
  data$type <- 23L
  data$filter <- filter
  data$correlations <- data$filter$acf
  data$argumentsList <- list(correlations = data$correlations, filterLength = filter$len)
  data$possibleLengths <- (data$filter$len + 2L):data$n
  data$defaultLengths <- (data$filter$len + 2L):data$n
  data$defaultIntervalSystem <- "dyaLen"
  data$defaultPenalty <- "weights"
  data$penaltyShift <- -data$filter$len
  
  data$ma <- .computeMA(data$filter$acf)
  data$key <- digest::sha1(list("hjsmurfLR", data$ma), digits = 6)
  data$rand.gen <- .rand.genMDependent
  data$addData <- function(rand.gen, data) {
    data$y <- rand.gen(data)
    
    if (!is.numeric(data$y) || length(data$y) != data$n || any(!is.finite(data$y))) {
      stop("for parametric family 'hjsmurfLR' rand.gen must return a finite numeric vector ",
           "of length equal to the number of observations")
    }
    data
  }
  data$generateSignal <- function(data, intervalSystem) {
    0
  }
  data$MC <- function(data, n) {
    data$n <- n
    data$nq <- n
    data$possibleLengths <- (data$filter$len + 2L):data$n
    data$defaultLengths <- (data$filter$len + 2L):data$n
    data
  }
  data$save <- TRUE
  
  data
}

.family2Param <- function(data, y, filter = NULL, fit = NULL, startTime = 0,
                          thresholdLongSegment = 25L, localValue = stats::median,
                          localVar = function(data) sdrobnorm(data, lag = filter$len + 1L)^2,
                          regularization = 1, correlations, suppressWarningNoDeconvolution = FALSE, localList = NULL) {
  if (!is.null(y)) {
    if (!is.double(y) || any(!is.finite(y))) {
      stop("data must be a finite numerical vector")
    }
    
    if (length(y) != data$n) {
      stop("data must be of length n")
    }
    
    data$y <- y
  }
  
  if (is.null(filter)) {
    stop("filter must be given for family '2Param'")
  } 
  
  if (!is(filter, "lowpassFilter")) {
    stop("filter must be an object of class 'lowpassFilter'")
  }
  
  if (!is.numeric(startTime) || length(startTime) != 1 || !is.finite(startTime)) {
    stop("startTime must be a single finite numeric")
  }
  time <- startTime + seq(along = data$y) / filter$sr
  
  if (is.null(data$y)) {
    fit <- stepblock(value = 0, leftEnd = startTime, rightEnd = startTime + data$n / filter$sr, x0 = startTime)
    fit$var <- 1
  } else {
    if (is.null(fit)) {
      fit <- stepblock(value = 0, leftEnd = startTime, rightEnd = startTime + data$n / filter$sr,  x0 = startTime)
      fit$var <- 1
    } else {
      if (!is(fit, "stepblock")) {
        if (is.list(fit) && !is.null(fit$fit) && is(fit$fit, "stepblock")) {
          fit <- fit$fit
        } else {
          stop("argument 'fit' must be an object of class 'stepblock'")
        }
      } else {
        if (fit$leftEnd[1] > startTime || fit$rightEnd[length(fit$rightEnd)] < time[length(time)]) {
          stop("fit does not provide values for all time points;",
               " please verify that the correct startTime and fit are given")
        }
      }
    }
  }
  
  if (!is.numeric(thresholdLongSegment) || length(thresholdLongSegment) != 1 ||
      !is.finite(thresholdLongSegment)) {
    stop("thresholdLongSegment must be a single positive integer")
  }
  
  if (!is.integer(thresholdLongSegment)) {
    thresholdLongSegment <- as.integer(thresholdLongSegment + 1e-6)
  }
  
  if (thresholdLongSegment <= 0L) {
    stop("thresholdLongSegment must be a single positive integer")
  }
  
  if (!is.function(localValue) || length(names(formals(localValue))) == 0) {
    stop("localValue must be a function with at least one argument")
  }
  
  if (!is.function(localVar) || length(names(formals(localVar))) == 0) {
    stop("localVar must be a function with at least one argument")
  }
  
  if (!is.numeric(regularization) || any(!is.finite(regularization))) {
    stop("all entries of 'regularization' must be finite numerics")
  }
  
  if (missing(correlations)) {
    regu <- regularization[1:(filter$len + 1)]
    regu[is.na(regu)] <- 0
    correlations <- filter$acf + regu
  } else {
    if (!is.numeric(correlations) || !all(is.finite(correlations))) {
      stop("correlations must be a finite numeric vector")
    }
    
    if (correlations[length(correlations)] == 0) {
      stop("the last element of correlations should not be zero")
    }
  }
  
  if (!is.logical(suppressWarningNoDeconvolution) || length(suppressWarningNoDeconvolution) != 1 || 
      is.na(suppressWarningNoDeconvolution)) {
    stop("suppressWarningNoDeconvolution must be a single logical (not NA)")
  }
  
  data$save <- identical(thresholdLongSegment, 25L) && identical(localValue, stats::median) &&
    identical(localVar, function(data) sdrobnorm(data, lag = filter$len + 1L)^2, ignore.environment = TRUE)
  
  if (is.null(data$y)) {
    add <- integer(0L)
  } else {
    # find long segments and reestimate their values
    shiftStart <- filter$len / filter$sr
    shiftEnd <- filter$len / filter$sr
    tolerance <- 1e-6 / filter$sr
    
    value <- numeric(length(fit$value))
    var <- numeric(length(fit$value))
    
    for (i in seq(along = fit$leftEnd)) {
      indices <- .whichIndices(time, fit$leftEnd[i] + shiftStart - tolerance,
                               fit$rightEnd[i] - shiftEnd + tolerance)
      
      if (length(indices) >= thresholdLongSegment) {
        est <- localValue(data$y[indices])
        if (!is.numeric(est) || length(est) != 1 || !is.finite(est)) {
          stop("localValue must return a single finite numeric")
        }
        value[i] <- est
        
        est <- localVar(data$y[indices])
        if (!is.numeric(est) || length(est) != 1 || !is.finite(est) || est <= 0) {
          stop("localVar must return a single finite positive numeric")
        }
        var[i] <- est
      } else {
        value[i] <- NA
        var[i] <- NA
      }
    }
    
    gridSize <- 1 / filter$sr
    
    jumps <- integer(0L)
    jumpsD <- integer(0L)
    jumps2 <- integer(0L)
    add <- integer(0L)
    addLeft <- integer(0L)
    addRight <- integer(0L)
    noDeconvolution <- logical(0L)
    
    if (length(value) == 1L && is.na(value)) {
      add <- 1:data$n
      addLeft <- as.integer((fit$leftEnd[1] - startTime) * filter$sr + 1e-6)
      addRight <- as.integer((fit$rightEnd[1] - startTime) * filter$sr + 1e-6)
      noDeconvolution <- TRUE
    }
    
    i <- 1L
    while (i < length(value)) {
      if (!is.na(value[i]) && !is.na(value[i + 1])) {
        # jump
        indices <- .whichIndices(time, fit$leftEnd[i + 1] - shiftEnd + tolerance,
                                 fit$rightEnd[i] + shiftStart - tolerance)
        gridIndices <- .whichIndices(time, fit$leftEnd[i + 1] - shiftEnd - tolerance, 
                                     fit$rightEnd[i] + tolerance)
        cp <- lowpassFilter::.deconvolveJump(seq(time[gridIndices[1]], time[gridIndices[length(gridIndices)]],
                                                 1 / filter$sr),
                                             data$y[indices], time[indices],
                                             as.numeric(value[i]), as.numeric(value[i + 1]),
                                             filter$number, filter$list, correlations)
        
        jumps[length(jumps) + 1] <- cp
        jumpsD[length(jumpsD) + 1] <- fit$rightEnd[i]
      } else if (!is.na(value[i]) && is.na(value[i + 1]) && i + 1 != length(value) && !is.na(value[i + 2])) {
        # isolated peak
        newAdd <- .whichIndices(time, fit$leftEnd[i + 1] - shiftEnd - tolerance,
                                fit$rightEnd[i + 1] + tolerance)
        add[length(add) + seq(along = newAdd)] <- newAdd
        addLeft[length(addLeft) + 1] <- as.integer((fit$leftEnd[i + 1] - startTime) * filter$sr + 1e-6)
        addRight[length(addRight) + 1] <- as.integer((fit$rightEnd[i + 1] - startTime) * filter$sr + 1e-6)
        jumps2[length(jumps2) + 1:2] <- c(fit$rightEnd[i], fit$rightEnd[i + 1L])
        noDeconvolution[length(noDeconvolution) + 1] <- FALSE
        i <- i + 1L
      } else {
        # deconvolution impossible
        if (!suppressWarningNoDeconvolution) {
          suppressWarningNoDeconvolution <- TRUE
          warning("at least one segment could not be deconvoled ",
                  "since two successive short segments (or a short segment at the begin or end) occured")
        }
        
        if (is.na(value[i])) {
          start <- fit$leftEnd[i] - shiftEnd - tolerance
          startIndex <- i
          if (i != 1L) {
            stop("unexpected result")
          }
        } else {
          start <- fit$leftEnd[i + 1] - shiftEnd - tolerance
          startIndex <- i + 1L
        }
        
        jumps2[length(jumps2) + 1] <- fit$rightEnd[i]
        
        while (is.na(value[i + 1])) {
          i <- i + 1L
          if (i == length(value)) {
            break
          }
          jumps2[length(jumps2) + 1] <- fit$rightEnd[i]
        }
        
        newAdd <- .whichIndices(time, start, fit$rightEnd[i] + tolerance)
        add[length(add) + seq(along = newAdd)] <- newAdd
        addLeft[length(addLeft) + 1] <- as.integer((fit$leftEnd[startIndex] - startTime) * filter$sr + 1e-6)
        addRight[length(addRight) + 1] <- as.integer((fit$rightEnd[i] - startTime) * filter$sr + 1e-6)
        noDeconvolution[length(noDeconvolution) + 1] <- TRUE
      }
      i <- i + 1L
    }
    
    fit <- stepblock(value = value, leftEnd = c(startTime, sort(c(jumps, jumps2))),
                     rightEnd = c(sort(c(jumps, jumps2)), time[length(time)]), x0 = startTime)
    fit$var <- var
  }
  
  data$type <- 100L
  data$fit <- fit
  data$filter <- filter
  data$thresholdLongSegment <- thresholdLongSegment
  data$localValue <- localValue
  data$localVar <- localVar
  data$startTime <- startTime
  data$add <- add
  data$localList <- localList
  
  if (is.null(data$y)) {
    data$argumentsList <- list()
  } else {
    data$time <- time
    data$addLeft <- addLeft
    data$addRight <- addRight
    data$jumps <- jumps
    data$jumpsD <- jumpsD
    data$noDeconvolution <- noDeconvolution
    data$tolerance <- tolerance
    
    val <- numeric(length(data$y))
    var <- numeric(length(data$y))
    stepFilter <- data$filter$truncatedStepfun(1:data$filter$len / data$filter$sr)
    stepAcf <- data$filter$acAntiderivative(1:data$filter$len / data$filter$sr, 0)
    left <- as.integer((data$fit$leftEnd - data$startTime) * data$filter$sr + 1e-6)
    right <- as.integer((data$fit$rightEnd - data$startTime) * data$filter$sr + 1e-6)
    
    val[left[1]:right[1]] <- data$fit$value[1]
    var[left[1]:right[1]] <- data$fit$var[1]
    
    for (i in seq(along = left)[-1]) {
      end <- min(left[i] + data$filter$len, right[i])
      number <- min(end - left[i], length(stepFilter))
      val[(left[i] + 1):end] <- data$fit$value[i - 1] * (1 - stepFilter[1:number]) +
        data$fit$value[i] * stepFilter[1:number]
      var[(left[i] + 1):end] <- data$fit$var[i - 1] * (1 - stepAcf[1:number]) + data$fit$var[i] * stepAcf[1:number]
      if ((left[i] + data$filter$len) <= right[i]) {
        val[(left[i] + data$filter$len):right[i]] <- data$fit$value[i]
        var[(left[i] + data$filter$len):right[i]] <- data$fit$var[i]
      }
    }
    val[add] <- NA
    var[add] <- NA
    
    Tobs <- (data$y - val)^2
    data$argumentsList <- list(obs = data$y, T0 = var, Tobs = Tobs, value = data$fit$value, var = data$fit$var,
                               filterLength = data$filter$len)
  }
  
  data$possibleLengths <- 1:data$n
  data$defaultLengths <- 1:min(data$n, 65L)
  data$defaultIntervalSystem <- "all"
  data$defaultPenalty <- "weights"
  data$penaltyShift <- 0
  data$ma <- .computeMA(data$filter$acf)
  
  # unname to be compatibel to old versions of digest::sha1
  data$key <- digest::sha1(list("2Param", filter$type, unname(filter$param),
                                filter$sr, filter$len), digits = 6)

  data$rand.gen <- .rand.genMDependent
  data$addData <- function(rand.gen, data) {
    data$y <- rand.gen(data)
    
    if (!is.numeric(data$y) || length(data$y) != data$n || any(!is.finite(data$y))) {
      stop("for parametric family '2Param' rand.gen must return a finite numeric vector ",
           "of length equal to the number of observations")
    }
    
    indices <- data$filter$len:(data$n - data$filter$len)
    if (length(indices) >= data$thresholdLongSegment) {
      est <- data$localValue(data$y[indices])
      if (!is.numeric(est) || length(est) != 1 || !is.finite(est)) {
        stop("localValue must return a single finite numeric")
      }
      
      var <- data$localVar(data$y[indices])
      if (!is.numeric(var) || length(var) != 1 || !is.finite(var) || var <= 0) {
        stop("localVar must return a single finite positive numeric")
      }
      data$y <- (data$y - est) / sqrt(var)
    } else {
      data$fit$value <- NA
    }
    
    # update data$argumentsList
    val <- rep(data$fit$value[1], length(data$y))
    var <- rep(data$fit$var[1], length(data$y))
    
    Tobs <- (data$y - val)^2
    data$argumentsList <- list(obs = data$y, T0 = var, Tobs = Tobs, value = data$fit$value, var = data$fit$var,
                               filterLength = data$filter$len)
    
    data
  }
  
  data$generateSignal <- function(data, intervalSystem) {
    if (intervalSystem$intervalSystem == "dyaPar") {
      stop("intervalSystem 'dyaPar' is not supported for this parametric family")
    }
    
    lengths <- intervalSystem$lengths
    left <- as.integer((data$fit$leftEnd - data$startTime) * data$filter$sr + 1e-6)
    right <- as.integer((data$fit$rightEnd - data$startTime) * data$filter$sr + 1e-6)
    
    if (is.null(data$localList)) {
      localList <- list()
      for (indexLen in seq(along = lengths)) {
        len <- lengths[indexLen]
        time <- 1:(len + data$filter$len - 1) / data$filter$sr
        cpLeft <- 0
        cpRight <- len / data$filter$sr
        
        Fleft <- data$filter$truncatedStepfun(time - cpLeft)
        Fright <- data$filter$truncatedStepfun(time - cpRight)
        v <- Fleft - Fright
        sumv2 <- sum(v^2)
        
        Fleft <- outer(time, time, function(i, j) data$filter$acAntiderivative(pmin(i, j) - cpLeft, abs(j - i)))
        Fright <- outer(time, time, function(i, j) data$filter$acAntiderivative(pmin(i, j) - cpRight, abs(j - i)))
        cor <- outer(time, time, function(i, j) data$filter$acfun(abs(j - i)))
        w <- Fleft - Fright
        sigmaL <- (cor - Fleft)
        sigmaR <- Fright
        vv <- outer(seq(along = time), seq(along = time), function(i, j) v[i] * v[j] / sum(v^2))
        diagW <- diag(w)
        matrixDiagW <- matrix(rep(diagW, length(diagW)), length(diagW))
        AL <- sum(diag(sigmaL) * diagW) - sum(vv * sigmaL * matrixDiagW)
        AR <- sum(diag(sigmaR) * diagW) - sum(vv * sigmaR * matrixDiagW)
        B <- sum(diagW^2) - sum(vv * w * matrixDiagW)
        
        w <- diagW
        sigmaL <- diag(sigmaL)
        sigmaR <- diag(sigmaR)
        
        Fleft <- 1 - data$filter$truncatedStepfun(time - cpLeft)
        Fright <- data$filter$truncatedStepfun(time - cpRight)
        
        localList[[indexLen]] = list(len = len, Fleft = Fleft, Fright = Fright, v = v, sumv2 = sumv2,
                                     sumSigmaL = AL, sumSigmaR = AR, sumW = B, w = w,
                                     sigmaL = sigmaL, sigmaR = sigmaR)
      }
    } else {
      if (is(data$localList, "localList") && identical(attr(data$localList, "method"), "2Param") && 
          identical(attr(data$localList, "filter"), data$filter, ignore.environment = TRUE) &&
          identical(attr(data$localList, "lengths"), lengths)) {
        localList <- data$localList
      } else {
        if (!is(data$localList, "localList")) {
          stop("localList is not an object of class 'localList'")
        }
        
        if (!identical(attr(data$localList, "method"), "2Param")) {
          stop("localList was created for a different method")
        }
        
        if (!identical(attr(data$localList, "filter"), data$filter, ignore.environment = TRUE)) {
          stop("localList was created for a different filter")
        }
        
        if (!identical(attr(data$localList, "lengths"), lengths)) {
          stop("localList was created for different lengths")
        }
      }
    }
    
    list(n = data$n, lengths = lengths, left = left, right = right - 1L, filterLength = data$filter$len,
         argumentsListLocal = localList)
  }

  data$MC <- function(data, n) {
    data$n <- n
    data$nq <- n
    data$possibleLengths <- 1:data$n
    data$defaultLengths <- 1:min(data$n, 65L)
    data$fit <- stepblock(value = 0, leftEnd = startTime,
                          rightEnd = startTime + data$n / data$filter$sr, x0 = startTime)
    data$fit$var <- 1
    data$localList <- NULL
    data
  }
  
  data
}

.familyLR <- function(data, y, filter = NULL, fit = NULL, startTime = 0,
                      thresholdLongSegment = 10L, localValue = stats::median,
                      regularization = 1, correlations, suppressWarningNoDeconvolution = FALSE, localList = NULL) {
  if (!is.null(y)) {
    if (!is.double(y) || any(!is.finite(y))) {
      stop("data must be a finite numerical vector")
    }
    
    if (length(y) != data$n) {
      stop("data must be of length n")
    }
    
    data$y <- y
  }
  
  if (is.null(filter)) {
    stop("filter must be given for family 'LR'")
  } 
  
  if (!is(filter, "lowpassFilter")) {
    stop("filter must be an object of class 'lowpassFilter'")
  }
  
  if (!is.numeric(startTime) || length(startTime) != 1 || !is.finite(startTime)) {
    stop("startTime must be a single finite numeric")
  }
  time <- startTime + seq(along = data$y) / filter$sr
  
  if (is.null(data$y)) {
    fit <- stepblock(value = 0, leftEnd = startTime, rightEnd = startTime + data$n / filter$sr, x0 = startTime)
  } else {
    if (is.null(fit)) {
      fit <- stepblock(value = 0, leftEnd = startTime, rightEnd = startTime + data$n / filter$sr, x0 = startTime)
    } else {
      if (!is(fit, "stepblock")) {
        if (is.list(fit) && !is.null(fit$fit) && is(fit$fit, "stepblock")) {
          fit <- fit$fit
        } else {
          stop("argument 'fit' must be an object of class 'stepblock'")
        }
      } else {
        if (fit$leftEnd[1] > startTime || fit$rightEnd[length(fit$rightEnd)] < time[length(time)]) {
          stop("fit does not provide values for all time points;",
               " please verify that the correct startTime and fit are given")
        }
      }
    }
  }
  
  if (!is.numeric(thresholdLongSegment) || length(thresholdLongSegment) != 1 ||
      !is.finite(thresholdLongSegment)) {
    stop("thresholdLongSegment must be a single positive integer")
  }
  
  if (!is.integer(thresholdLongSegment)) {
    thresholdLongSegment <- as.integer(thresholdLongSegment + 1e-6)
  }
  
  if (thresholdLongSegment <= 0L) {
    stop("thresholdLongSegment must be a single positive integer")
  }
  
  if (!is.function(localValue) || length(names(formals(localValue))) == 0) {
    stop("localValue must be a function with at least one argument")
  }
  
  if (!is.numeric(regularization) || any(!is.finite(regularization))) {
    stop("all entries of 'regularization' must be finite numerics")
  }
  
  if (missing(correlations)) {
    regu <- regularization[1:(filter$len + 1)]
    regu[is.na(regu)] <- 0
    correlations <- filter$acf + regu
  } else {
    if (!is.numeric(correlations) || !all(is.finite(correlations))) {
      stop("correlations must be a finite numeric vector")
    }
    
    if (correlations[length(correlations)] == 0) {
      stop("the last element of correlations should not be zero")
    }
  }
  
  if (!is.logical(suppressWarningNoDeconvolution) || length(suppressWarningNoDeconvolution) != 1 || 
      is.na(suppressWarningNoDeconvolution)) {
    stop("suppressWarningNoDeconvolution must be a single logical (not NA)")
  }
  
  data$save <- identical(thresholdLongSegment, 10L) && identical(localValue, stats::median)
  
  if (is.null(data$y)) {
    add <- integer(0L)
  } else {
    # find long segments and reestimate their values
    shiftStart <- filter$len / filter$sr
    shiftEnd <- filter$len / filter$sr
    tolerance <- 1e-6 / filter$sr
    
    value <- numeric(length(fit$value))
    
    for (i in seq(along = fit$leftEnd)) {
      indices <- .whichIndices(time, fit$leftEnd[i] + shiftStart - tolerance,
                                              fit$rightEnd[i] - shiftEnd + tolerance)
      
      if (length(indices) >= thresholdLongSegment) {
        est <- localValue(data$y[indices])
        if (!is.numeric(est) || length(est) != 1 || !is.finite(est)) {
          stop("localValue must return a single finite numeric")
        }
        value[i] <- est
      } else {
        value[i] <- NA
      }
    }
    
    gridSize <- 1 / filter$sr
    
    jumps <- integer(0L)
    jumpsD <- integer(0L)
    jumps2 <- integer(0L)
    add <- integer(0L)
    addLeft <- integer(0L)
    addRight <- integer(0L)
    noDeconvolution <- logical(0L)
    
    if (length(value) == 1L && is.na(value)) {
      add <- 1:data$n
      addLeft <- as.integer((fit$leftEnd[1] - startTime) * filter$sr + 1e-6)
      addRight <- as.integer((fit$rightEnd[1] - startTime) * filter$sr + 1e-6)
      noDeconvolution <- TRUE
    }
    
    i <- 1L
    while (i < length(value)) {
      if (!is.na(value[i]) && !is.na(value[i + 1])) {
        # jump
        indices <- .whichIndices(time, fit$leftEnd[i + 1] - shiftEnd + tolerance,
                                 fit$rightEnd[i] + shiftStart - tolerance)
        gridIndices <- .whichIndices(time, fit$leftEnd[i + 1] - shiftEnd - tolerance, 
                                     fit$rightEnd[i] + tolerance)
        cp <- lowpassFilter::.deconvolveJump(seq(time[gridIndices[1]], time[gridIndices[length(gridIndices)]],
                                                 1 / filter$sr),
                              data$y[indices], time[indices],
                              as.numeric(value[i]), as.numeric(value[i + 1]),
                              filter$number, filter$list, correlations)
        
        jumps[length(jumps) + 1] <- cp
        jumpsD[length(jumpsD) + 1] <- fit$rightEnd[i]
      } else if (!is.na(value[i]) && is.na(value[i + 1]) && i + 1 != length(value) && !is.na(value[i + 2])) {
        # isolated peak
        newAdd <- .whichIndices(time, fit$leftEnd[i + 1] - shiftEnd - tolerance,
                                fit$rightEnd[i + 1] + tolerance)
        add[length(add) + seq(along = newAdd)] <- newAdd
        addLeft[length(addLeft) + 1] <- as.integer((fit$leftEnd[i + 1] - startTime) * filter$sr + 1e-6)
        addRight[length(addRight) + 1] <- as.integer((fit$rightEnd[i + 1] - startTime) * filter$sr + 1e-6)
        jumps2[length(jumps2) + 1:2] <- c(fit$rightEnd[i], fit$rightEnd[i + 1L])
        noDeconvolution[length(noDeconvolution) + 1] <- FALSE
        i <- i + 1L
      } else {
        # deconvolution impossible
        if (!suppressWarningNoDeconvolution) {
          suppressWarningNoDeconvolution <- TRUE
          warning("at least one segment could not be deconvoled ",
                  "since two successive short segments (or a short segment at the begin or end) occured")
        }
        
        if (is.na(value[i])) {
          start <- fit$leftEnd[i] - shiftEnd - tolerance
          startIndex <- i
          if (i != 1L) {
            stop("unexpected result")
          }
        } else {
          start <- fit$leftEnd[i + 1] - shiftEnd - tolerance
          startIndex <- i + 1L
        }
        
        jumps2[length(jumps2) + 1] <- fit$rightEnd[i]
        
        while (is.na(value[i + 1])) {
          i <- i + 1L
          if (i == length(value)) {
            break
          }
          jumps2[length(jumps2) + 1] <- fit$rightEnd[i]
        }
        
        newAdd <- .whichIndices(time, start, fit$rightEnd[i] + tolerance)
        add[length(add) + seq(along = newAdd)] <- newAdd
        addLeft[length(addLeft) + 1] <- as.integer((fit$leftEnd[startIndex] - startTime) * filter$sr + 1e-6)
        addRight[length(addRight) + 1] <- as.integer((fit$rightEnd[i] - startTime) * filter$sr + 1e-6)
        noDeconvolution[length(noDeconvolution) + 1] <- TRUE
      }
      i <- i + 1L
    }
    
    fit <- stepblock(value = value, leftEnd = c(startTime, sort(c(jumps, jumps2))),
                     rightEnd = c(sort(c(jumps, jumps2)), time[length(time)]), x0 = startTime)
  }
  
  data$type <- 102L
  data$fit <- fit
  data$filter <- filter
  data$thresholdLongSegment <- thresholdLongSegment
  data$localValue <- localValue
  data$startTime <- startTime
  data$add <- add
  data$localList <- localList
  
  if (is.null(data$y)) {
    data$argumentsList <- list()
  } else {
    data$time <- time
    data$addLeft <- addLeft
    data$addRight <- addRight
    data$jumps <- jumps
    data$jumpsD <- jumpsD
    data$noDeconvolution <- noDeconvolution
    data$tolerance <- tolerance
    
    covariances <- data$filter$acf * sdrobnorm(data$y, lag = data$filter$len + 1L)^2
    covariances[1] <- covariances[1] * 2
    
    val <- numeric(length(data$y))
    stepFilter <- data$filter$truncatedStepfun(1:data$filter$len / data$filter$sr)
    left <- as.integer((data$fit$leftEnd - data$startTime) * data$filter$sr + 1e-6)
    right <- as.integer((data$fit$rightEnd - data$startTime) * data$filter$sr + 1e-6)
    
    val[left[1]:right[1]] <- data$fit$value[1]
    for (i in seq(along = left)[-1]) {
      end <- min(left[i] + data$filter$len, right[i])
      number <- min(end - left[i], length(stepFilter))
      val[(left[i] + 1):end] <- data$fit$value[i - 1] * (1 - stepFilter[1:number]) +
        data$fit$value[i] * stepFilter[1:number]
      if ((left[i] + filter$len) <= right[i]) {
        val[(left[i] + filter$len):right[i]] <- data$fit$value[i]
      }
    }
    val[data$add] <- NA
    
    obs0 <- data$y - val
    data$argumentsList <- list(obs = data$y, obs0 = obs0, value = data$fit$value, covariances = covariances,
                               filterLength = data$filter$len)
  }
  
  data$possibleLengths <- 1:data$n
  data$defaultLengths <- 1:min(data$n, 20L)
  data$defaultIntervalSystem <- "all"
  data$defaultPenalty <- "weights"
  data$penaltyShift <- 0
  data$ma <- .computeMA(data$filter$acf)
  
  # unname to be compatibel to old versions of digest::sha1
  data$key <- digest::sha1(list("LR", filter$type, unname(filter$param),
                                filter$sr, filter$len), digits = 6)
  
  data$rand.gen <- .rand.genMDependent
  data$addData <- function(rand.gen, data) {
    data$y <- rand.gen(data)
    
    if (!is.numeric(data$y) || length(data$y) != data$n || any(!is.finite(data$y))) {
      stop("for parametric family 'LR' rand.gen must return a finite numeric vector ",
           "of length equal to the number of observations")
    }
    
    indices <- data$filter$len:(data$n - data$filter$len)
    if (length(indices) >= data$thresholdLongSegment) {
      est <- data$localValue(data$y[indices])
      if (!is.numeric(est) || length(est) != 1 || !is.finite(est)) {
        stop("localValue must return a single finite numeric")
      }

      data$y <- data$y - est
    } else {
      data$fit$value <- NA
    }
    
    # update data$argumentsList
    covariances <- data$filter$acf * sdrobnorm(data$y, lag = data$filter$len + 1L)^2
    covariances[1] <- covariances[1] * 2
    
    data$argumentsList <- list(obs = data$y, obs0 = data$y, value = data$fit$value,
                               covariances = covariances, filterLength = data$filter$len)
    
    data
  }
  
  data$generateSignal <- function(data, intervalSystem) {
    if (intervalSystem$intervalSystem == "dyaPar") {
      stop("intervalSystem 'dyaPar' is not supported for this parametric family")
    }
    
    lengths <- intervalSystem$lengths
    left <- as.integer((data$fit$leftEnd - data$startTime) * data$filter$sr + 1e-6)
    right <- as.integer((data$fit$rightEnd - data$startTime) * data$filter$sr + 1e-6)
    
    if (is.null(data$localList)) {
      correlations <- data$filter$acf
      correlations[1] <- correlations[1] + 1
      
      localList <- list()
      for (indexLen in seq(along = lengths)) {
        len <- lengths[indexLen]
        time <- 1:(len + data$filter$len - 1) / data$filter$sr
        cpLeft <- 0
        cpRight <- len / data$filter$sr
        
        m <- min(len + data$filter$len - 1, length(correlations) - 1L)
        
        A <- matrix(0, len + data$filter$len - 1, len + data$filter$len - 1)
        for (i in 1:(len + data$filter$len - 2)) {
          A[i, i] <- correlations[1]
          A[i, i + 1:min(m, len + data$filter$len - 1 - i)] <- 
            correlations[2:min(m + 1, len + data$filter$len - 1 - i + 1)]
          A[i + 1:min(m, len + data$filter$len - 1 - i), i] <- 
            correlations[2:min(m + 1, len + data$filter$len - 1 - i + 1)]
        }
        A[len + data$filter$len - 1, len + data$filter$len - 1] <- correlations[1]
        
        Fleft <- data$filter$truncatedStepfun(time - cpLeft)
        Fright <- data$filter$truncatedStepfun(time - cpRight)
        v <- Fleft - Fright
        sol <- solve(A, v)
        vtAv <- sum(v * sol)
        
        Fleft <- 1 - Fleft
        
        localList[[indexLen]] = list(len = len, Fleft = Fleft, Fright = Fright, v = v, sol = sol, vtAv = vtAv)
      }
    } else {
      if (is(data$localList, "localList") && identical(attr(data$localList, "method"), "LR") && 
          identical(attr(data$localList, "filter"), data$filter, ignore.environment = TRUE) &&
          identical(attr(data$localList, "lengths"), lengths)) {
        localList <- data$localList
      } else {
        if (!is(data$localList, "localList")) {
          stop("localList is not an object of class 'localList'")
        }
        
        if (!identical(attr(data$localList, "method"), "LR")) {
          stop("localList was created for a different method")
        }
        
        if (!identical(attr(data$localList, "filter"), data$filter, ignore.environment = TRUE)) {
          stop("localList was created for a different filter")
        }
        
        if (!identical(attr(data$localList, "lengths"), lengths)) {
          stop("localList was created for different lengths")
        }
      }
    }
    
    list(n = data$n, lengths = lengths, left = left, right = right - 1L, filterLength = data$filter$len,
         argumentsListLocal = localList)
  }
  
  data$MC <- function(data, n) {
    data$n <- n
    data$nq <- n
    data$possibleLengths <- 1:data$n
    data$defaultLengths <- 1:min(data$n, 20L)
    data$fit <- stepblock(value = 0, leftEnd = startTime, 
                          rightEnd = startTime + data$n / data$filter$sr, x0 = startTime)
    data$localList <- NULL
    data
  }
  
  data
}

.rand.genMDependent <- function(data) {
  kern <- c(1, data$ma)
  z <- rnorm(data$n + length(kern) - 1, sd = 1)
  lowpassFilter::.convolve(z, kern) / sqrt(sum(kern^2))
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

# indices such that time > start and time < end
# which(time > start & time < end)
.whichIndices <- function(time, start, end) {
  if (start > end) {
    return(integer(0))
  }
  
  left <- 1L
  right <- length(time)
  m <- as.integer((left + right) / 2 + 0.1)
  
  while (left != right) {
    if (time[m] > start) {
      right <- m
    } else {
      left <- m + 1
    }
    m <- as.integer((left + right) / 2 + 0.1)
  }
  
  start <- left
  
  right <- length(time)
  m <- as.integer((left + right) / 2 + 0.6)
  
  while (left != right) {
    if (time[m] < end) {
      left <- m
    } else {
      right <- m - 1
    }
    m <- as.integer((left + right) / 2 + 0.6)
  }
  
  end <- left
  if (start <= end) {
    indices <- start:end
  } else {
    indices <- integer(0)
  }
  indices
}
