monteCarloSimulation <- function(n, r = 1e4L, family = NULL, intervalSystem = NULL, lengths = NULL,
                                 penalty = NULL, output = c("vector", "maximum"), seed = n,
                                 rand.gen = NULL, messages = NULL, ...) {
  if (missing(n)) {
    stop("number of observations 'n' must be given to perform Monte-Carlo simulations")
  }
  
  data <- .parametricFamily(family = family, n = n, nq = n, ...)
  data <- data$MC(data, data$n)
  
  output <- match.arg(output)
  
  if (output == "vector") {
    if (!is.null(lengths)) {
      warning("argument 'lengths' is given, but will be ignored for output 'vector'")
    }
    
    if (!is.null(penalty) && penalty != "none") {
      warning("argument 'penalty' is given, but will be ignored for output 'vector'")
    }
    
    intervalSystem <- .intervalSystem(intervalSystem = intervalSystem, lengths = NULL, data = data)
    penalty <- "none"
  } else {
    intervalSystem <- .intervalSystem(intervalSystem = intervalSystem, lengths = lengths, data = data)
    if (is.null(penalty)) {
      penalty <- data$defaultPenalty
    }
    if (penalty == "weights") {
      penalty <- "none"
    }
    penalty <- match.arg(penalty, c("sqrt", "log", "none"))
  }
  
  data$signal <- data$generateSignal(data = data, intervalSystem = intervalSystem)
  
  .monteCarloSimulation(data = data, intervalSystem = intervalSystem, penalty = penalty, r = r,
                        seed = seed, rand.gen = rand.gen, messages = messages, output = output)
}

.monteCarloSimulation <- function(data, intervalSystem, penalty, r = NULL, seed = data$n,
                                  rand.gen = NULL, messages = NULL, output) {
  if (!is.numeric(r) || length(r) != 1 || !is.finite(r)) {
    stop("r must be a single positive integer")
  }
  
  if (!is.integer(r)) {
    r <- as.integer(r + 1e-6)
  }
  
  if (r < 1L) {
    stop("r must be a single positive integer")
  }
  
  if (length(seed) == 1 && seed == "no") {
    # set no seed
  } else {
    set.seed(seed)
  }
  
  if (!is.null(messages)) {
    if (!is.function(messages)) {
      if (!is.numeric(messages) || length(messages) != 1 || !is.finite(messages)) {
        stop("messages must be a single positive integer")
      }
      
      if (!is.integer(messages)) {
        messages <- as.integer(messages + 1e-6)
      }
      
      if (messages < 1L) {
        stop("messages must be a single positive integer")
      }
      
      each <- messages
      messages <- function(i, r, each) {
        if (i %% each == 0) {
          message(paste(i, " of ", r, " simulations are completed.", sep = ""))
        }
        FALSE
      }
    } else {
      each <- NULL
    }
  } else {
    each <- NULL
    messages <- function(i, r, each) {FALSE}
  }
  
  if (output == "vector") {
    stat <- matrix(0, length(intervalSystem$lengths), r)
    keyList <- list(data$n, data$key, intervalSystem$intervalSystem)
    class(stat) <- c("MCSimulationVector", class(stat))
  } else {
    stat <- numeric(r)
    keyList <- list(data$n, data$key, intervalSystem$intervalSystem, intervalSystem$lengths, penalty)
    class(stat) <- c("MCSimulationMaximum", class(stat))
  }  
  
  if (is.null(rand.gen)) {
    rand.gen <- data$rand.gen
    save <- data$save
  } else {
    save <- FALSE
    if (!is.function(rand.gen) || !identical(names(formals(rand.gen)), "data")) {
      stop("rand.gen must be a function with a single argument data")
    }
  }
  
  for (i in 1:r) {
    if (messages(i, r, each)) {
      stop("User interrupt!")
    }
    
    data <- data$addData(rand.gen, data)

    if (output == "vector") {
      stat[, i] <- .computeStat(signal = data$signal, data = data, intervalSystem = intervalSystem, penalty = penalty,
                                output = output)
    } else {
      stat[i] <- .computeStat(signal = data$signal, data = data, intervalSystem = intervalSystem, penalty = penalty,
                              output = output)
    }
  }
  
  attr(stat, "keyList") <- keyList
  attr(stat, "key") <- digest::digest(keyList)
  attr(stat, "n") <- data$n
  attr(stat, "lengths") <- intervalSystem$lengths
  attr(stat, "save") <- save
  
  stat
}
