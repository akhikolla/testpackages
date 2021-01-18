computeStat  <- function(y, signal = 0, family = NULL, intervalSystem = NULL, lengths = NULL,
                         penalty = NULL, nq = length(y), output = c("list", "vector", "maximum"), ...) {
  if (missing(y)) {
    stop("argument 'y' must be given")
  }
  
  data <- .parametricFamily(family = family, y = y, n = length(y), nq = nq, ...)
  intervalSystem <- .intervalSystem(intervalSystem = intervalSystem, lengths = lengths, data = data)
  if (is.null(penalty)) {
    penalty <- data$defaultPenalty
  }
  penalty <- match.arg(penalty, c("sqrt", "log", "none", "weights"))
  if (penalty == "weights") {
    penalty <- "none"
  }
  output <- match.arg(output)
  
  if (data$type == 100L || data$type == 102L) {
    if (!identical(signal, 0)) {
      warning("argument 'signal' will be ignored for this parametric family,",
              " please use the parametric family argument 'fit' instead")
    }
    return(.computeStat(signal = data$generateSignal(data = data, intervalSystem = intervalSystem),
                        data = data, intervalSystem = intervalSystem, penalty = penalty, output = output)) 
  }

  if (is.list(signal)) {
    if (is.null(signal$leftIndex) || is.null(signal$rightIndex) || is.null(signal$value)) {
      stop("if signal is a list it must contain 'leftIndex', 'rightIndex' and 'value'")
    }
    
    if (!is.numeric(signal$leftIndex) || !all(is.finite(signal$leftIndex))) {
      stop("signal$leftIndex must be a integer vector")
    }
    
    if (!is.numeric(signal$rightIndex) || !all(is.finite(signal$rightIndex))) {
      stop("signal$rightIndex must be a integer vector")
    }
    
    if (!is.numeric(signal$value) || !all(is.finite(signal$value))) {
      stop("signal$value must be a finite numeric vector")
    }
    
    if (length(signal$value) != length(signal$leftIndex) ||
        length(signal$value) != length(signal$rightIndex)) {
      stop("signal$leftIndex, signal$rightIndex and signal$value must be of the same length")
    }
    
    if (!is.integer(signal$leftIndex)) {
      signal$leftIndex <- as.integer(signal$leftIndex + 1e-6)
    }
    
    if (!is.integer(signal$rightIndex)) {
      signal$rightIndex <- as.integer(signal$rightIndex + 1e-6)
    }
    
    if (data$n != signal$rightIndex[length(signal$rightIndex)] - signal$leftIndex[1] + 1L) {
      stop("length of observations must be the same as the length described by the signal")
    }
    
    if (signal$leftIndex[1] != 1L) {
      stop("signal$leftIndex must start at 1")
    }
  } else {
    if (!is.numeric(signal) || length(signal) != 1 || !is.finite(signal)) {
      stop("signal must be either a list or a single finite numeric")
    }
    
    if (signal == 0) {
      return(.computeStat(signal = signal, data = data, intervalSystem = intervalSystem, penalty = penalty,
                          output = output)) 
    } else {
      signal <- list(leftIndex = 1L, rightIndex = data$n, value = as.numeric(signal))
    }
  }
  
  signal <- list(leftIndex = signal$leftIndex - 1L, rightIndex = signal$rightIndex - 1L, # C style
                 value = as.numeric(signal$value))
  
  .computeStat(signal = signal, data = data, intervalSystem = intervalSystem, penalty = penalty,
               output = output)  
}

.computeStat <- function(signal, data, intervalSystem, penalty, output) {
  if (is.list(signal)) {
    stat <- .callRoutines(observations = data$y, routineType = 1L, argumentsListRoutine = signal, 
                          dataType = data$type, argumentsListData = data$argumentsList,
                          intervalSystemType = intervalSystem$type,
                          argumentsListIntervalSystem = intervalSystem$argumentsList)[intervalSystem$lengths]
  } else {
    stat <- .callRoutines(observations = data$y, routineType = 0L, argumentsListRoutine = NULL, 
                          dataType = data$type, argumentsListData = data$argumentsList,
                          intervalSystemType = intervalSystem$type,
                          argumentsListIntervalSystem = intervalSystem$argumentsList)[intervalSystem$lengths]
  }
  
  finiteLengths <- is.finite(stat)
  stat <- stat[finiteLengths]
  
  stat <- switch(penalty,
                 sqrt = sqrt(2 * stat) - 
                   sqrt(2 * log(exp(1) * data$nq / (intervalSystem$lengths[finiteLengths] + data$penaltyShift))),
                 log = stat - log(exp(1) * data$nq / (intervalSystem$lengths[finiteLengths] + data$penaltyShift)),
                 none = stat)
  
  if (output == "list") {
    ret <- list(maximum = max(stat), stat = rep(-Inf, length(intervalSystem$lengths)), 
                lengths = intervalSystem$lengths)
    ret$stat[finiteLengths] <- stat    
  } else if (output == "vector") {
    ret <- rep(-Inf, length(intervalSystem$lengths))
    ret[finiteLengths] <- stat
  } else {
    ret <- max(stat)
  }
  
  ret        
}
