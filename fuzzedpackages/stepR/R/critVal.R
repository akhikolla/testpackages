critVal <- function(n, q = NULL, alpha = NULL, nq = 2L^(as.integer(log2(n) + 1e-12) + 1L) - 1L, family = NULL,
                    intervalSystem = NULL, lengths = NULL, penalty = NULL, weights = NULL,
                    stat = NULL, r = 1e4, output = c("vector", "value"), options = NULL, ...) {
  if (missing(n)) {
    stop("number of observations 'n' must be given for computing critical values")
  }

  .RemoveAdditionalArgsPF <- function(family, n, nq, ..., seed, rand.gen, messages)
    .parametricFamily(family = family, n = n, nq = nq, ...)
  data <- .RemoveAdditionalArgsPF(family = family, n = n, nq = nq, ...)
  intervalSystem <- .intervalSystem(intervalSystem = intervalSystem, lengths = lengths, data = data)
  data$signal <- data$generateSignal(data = data, intervalSystem = intervalSystem)
  output <- match.arg(output)
  
  .RemoveAdditionalArgsCV <- function(q, data, output, intervalSystem, penalty, alpha, stat, r,
                                      weights, options, ..., sd, covariances, correlations, filter,
                                      fit, startTime, thresholdLongSegment, localValue, localVar,
                                      regularization, suppressWarningNoDeconvolution, localList)
    .critVal(q = q, data = data, output = output, intervalSystem = intervalSystem, penalty = penalty,
             alpha = alpha, stat = stat, r = r, weights = weights, options = options, ...)
  .RemoveAdditionalArgsCV(q = q, data = data, output = output, intervalSystem = intervalSystem, penalty = penalty,
                          alpha = alpha, stat = stat, r = r, weights = weights, options = options, ...)
}
  
.critVal <- function(q = NULL, output, data, intervalSystem, penalty = NULL,
                     alpha = 0.05, stat = NULL, r = 1e4, weights = NULL,
                     options = NULL, ...) {
  if (is.null(q)) {
    if (is.null(alpha)) {
      alpha <- 0.5
      warning("Either q or alpha must be given. ",
              "Please see the documentation for detailed information on this choise. ",
              "By default we continue with alpha = 0.5 which is often a reasonable choice.")
    }
    
    if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
      stop("alpha must be a probability, i.e., a single numeric between 0 and 1")
    }
    
    if (!is.null(stat) && !is.null(attr(stat, "n"))) {
      if (!is.numeric(attr(stat, "n")) || length(attr(stat, "n")) != 1 || !is.finite(attr(stat, "n"))) {
        stop("attribute n of stat must be a single integer larger than or equal to n")
      }
      
      if (!is.integer(attr(stat, "n"))) {
        attr(stat, "n") <- as.integer(attr(stat, "n") + 1e-6)
      }
      
      if (attr(stat, "n") < data$n) {
        stop("stat has been simulated for less than n = ", data$n, " observations (", attr(stat, "n"), "). ",
             "The Monte-Carlo simulation must be given for at least n observations")
      }
      
      data$nq <- attr(stat, "n")
    }
  } else {
    if (!is.numeric(q) || !all(is.finite(q))) {
      stop("q must be a finite numeric, either a single value or a vector")
    }
    
    if (!is.null(attr(q, "n"))) {
      if (!is.numeric(attr(q, "n")) || length(attr(q, "n")) != 1 || !is.finite(attr(q, "n"))) {
        stop("attribute n of q must be a single integer larger than or equal to n")
      }
      
      if (!is.integer(attr(q, "n"))) {
        attr(q, "n") <- as.integer(attr(q, "n") + 1e-6)
      }
      
      if (attr(q, "n") < data$n) {
        stop("q has been computed for less than n = ", data$n, " observations (", attr(q, "n"), "). ",
             "Must have been computed for at least n observations")
      }
      
      data$nq <- attr(q, "n")
    }
  }
  
  if (is.null(penalty)) {
    penalty <- data$defaultPenalty
  }
  penalty <- match.arg(penalty, c("sqrt", "log", "none", "weights"))
  
  if (is.null(options)) {
    if ((data$type == 100L || data$type == 102L) && !all(intervalSystem$lengths %in% data$defaultLengths)) {
      options <- list(simulation = "matrixIncreased",
                      load = list(), save = list(), envir = .GlobalEnv, dirs = "stepR")
    } else {
      options <- list(simulation = "matrixIncreased",
                      save = list(workspace = c("vector", "vectorIncreased"),
                                  fileSystem = c("matrix", "matrixIncreased"),
                                  RDSfile = NULL, variable = NULL),
                      load = list(workspace = c("vector", "vectorIncreased", "matrix", "matrixIncreased"),
                                  fileSystem = c("vector", "vectorIncreased", "matrix", "matrixIncreased"),
                                  package = TRUE, RDSfile = NULL),
                      envir = .GlobalEnv, dirs = "stepR")
    }
  } else {
    if (!is.list(options) || !all(names(options) %in% c("simulation", "save", "load", "envir", "dirs"))) {
      stop("options must be a list and only entries 'simulation', 'save', 'load', 'envir', 'dirs' are allowed")
    }
    
    options$simulation <- match.arg(options$simulation, c("matrixIncreased", "vector", "vectorIncreased", "matrix"))
    
    if (is.null(options$save)) {
      options$save <- list(workspace = c("vector", "vectorIncreased"),
                           fileSystem = c("matrix", "matrixIncreased"),
                           RDSfile = NULL, variable = NULL)
    } else {
      if (!is.list(options$save) ||
          !all(names(options$save) %in% c("workspace", "fileSystem", "RDSfile", "variable"))) {
        stop("options$save must be a list and ",
             "only entries 'workspace', 'fileSystem', 'RDSfile', 'variable' are allowed")
      }
      
      if (!is.null(options$save$workspace)) {
        for (i in seq(along = options$save$workspace)) {
          options$save$workspace[i] <- match.arg(options$save$workspace[i],
                                                 c("matrixIncreased", "vector", "vectorIncreased", "matrix"))
        }
        
        options$save$workspace <- match.arg(options$save$workspace,
                                            c("matrixIncreased", "vector", "vectorIncreased", "matrix"),
                                            several.ok = TRUE)
      }
      
      if (!is.null(options$save$fileSystem)) {
        for (i in seq(along = options$save$fileSystem)) {
          options$save$fileSystem[i] <- match.arg(options$save$fileSystem[i],
                                                  c("matrixIncreased", "vector", "vectorIncreased", "matrix"))
        }
        
        options$save$fileSystem <- match.arg(options$save$fileSystem,
                                             c("matrixIncreased", "vector", "vectorIncreased", "matrix"),
                                             several.ok = TRUE)
      }
      
      if (!is.null(options$save$RDSfile)) {
        if (length(options$save$RDSfile) == 1) {
          options$save$RDSfile <- rep(options$save$RDSfile, 2)
        }
        
        if (length(options$save$RDSfile) != 2 || any(!is.character(options$save$RDSfile))) {
          stop("options$save$RDSfile must have length one or two and contain paths to RDS file(s)")
        }
      }
      
      if (!is.null(options$save$variable)) {
        if (length(options$save$variable) == 1) {
          options$save$variable <- rep(options$save$variable, 2)
        }
        
        if (length(options$save$variable) != 2 || any(!is.character(options$save$variable))) {
          stop("options$save$variable must have length one or two and be strings that specify variables")
        }
      }
    }
    
    if (is.null(options$load)) {
      options$load <- list(workspace = c("vector", "vectorIncreased", "matrix", "matrixIncreased"),
                           fileSystem = c("vector", "vectorIncreased", "matrix", "matrixIncreased"),
                           package = TRUE, RDSfile = NULL)
    } else {
      if (!is.list(options$load) || !all(names(options$load) %in% c("workspace", "fileSystem",
                                                                    "package", "RDSfile"))) {
        stop("options$load must be a list and only entries 'workspace', 'fileSystem', 'RDSfile' are allowed")
      }
      
      if (!is.null(options$load$workspace)) {
        for (i in seq(along = options$load$workspace)) {
          options$load$workspace[i] <- match.arg(options$load$workspace[i],
                                                 c("matrixIncreased", "vector", "vectorIncreased", "matrix"))
        }
        
        options$load$workspace <- match.arg(options$load$workspace,
                                            c("matrixIncreased", "vector", "vectorIncreased", "matrix"),
                                            several.ok = TRUE)
      }
      
      if (!is.null(options$load$fileSystem)) {
        for (i in seq(along = options$load$fileSystem)) {
          options$load$fileSystem[i] <- match.arg(options$load$fileSystem[i],
                                                  c("matrixIncreased", "vector", "vectorIncreased", "matrix"))
        }
        
        options$load$fileSystem <- match.arg(options$load$fileSystem,
                                             c("matrixIncreased", "vector", "vectorIncreased", "matrix"),
                                             several.ok = TRUE)
      }
      
      if (!is.null(options$load$package)) {
        if (length(options$load$package) != 1 || !is.logical(options$load$package)) {
          stop("options$load$package must be a single logical")
        }
      }
      
      if (!is.null(options$load$RDSfile)) {
        if (length(options$load$RDSfile) != 1 || any(!is.character(options$load$RDSfile))) {
          stop("options$load$RDSfile must be a path to a RDS file")
        }
      }
    }
    
    if (is.null(options$envir)) {
      options$envir <- .GlobalEnv
    } else {
      if (!is.environment(options$envir)) {
        stop("options$envir must be an environment")
      }
    }
    
    if (is.null(options$dirs)) {
      options$dirs <- "stepR"
    } else {
      if (length(options$dirs) != 1 || !is.character(options$dirs)) {
        stop("options$dirs must be a string")
      }
      
      if (options$dirs == "") {
        options$dirs <- NULL
      }
    }
    
    if ((data$type == 100L || data$type == 102L) && !all(intervalSystem$lengths %in% data$defaultLengths)) {
      options$save$workspace <- NULL
      options$save$fileSystem <- NULL
      options$load$workspace <- NULL
      options$load$fileSystem <- NULL
      options$load$package <- FALSE
    }
  }
  
  if (penalty == "weights") {
    if (output == "value") {
      stop("output 'value' is not possible for penalty 'weights'")
    }
    
    if (is.null(q)) {
      if (!is.null(stat) && !is(stat, "MCSimulationVector")) {
        stop("stat must be an object of class 'MCSimulationVector' for penalty 'weights'")
      }
      
      stat <- .getStat(stat = stat, data = data, intervalSystem = intervalSystem, r = r, options = options, ...)
      q <- .critValWeights(stat = stat, alpha = alpha, weights = weights)
      attr(q, "n") <- attr(stat, "n")
    }
  } else {
    if (!is.null(weights)) {
      warning("weights are given, but penalty is not 'weights', weights will be ignored")
    }
    
    if (is.null(q)) {
      stat <- .getVectorStat(stat = stat, data = data, intervalSystem = intervalSystem,
                             penalty = penalty, r = r, options = options, ...)

      q <- as.numeric(stats::quantile(stat, 1 - alpha, type = 1))
      attr(q, "n") <- attr(stat, "n")
    }
  }
  
  if (output == "value") {
    if (length(q) != 1L) {
      stop("q must be a single finite numeric or NULL for output 'value'")
    }
    
    if (is.null(attr(q, "n"))) {
      attr(q, "n") <- data$n
    }
  } else {
    q <- .critValVector(q = q, data = data, intervalSystem = intervalSystem, penalty = penalty)
  }
  
  q
}

.getStat <- function(stat, data, intervalSystem, r = NULL, options = NULL, ...) {
  if (is.null(stat)) {
    stat <- .loadMatrix(data = data, intervalSystem = intervalSystem, r = r, options = options, ...)
    if (is.null(stat)) {
      stop("options$simulation must be 'matrix' or 'matrixIncreased' for penalty 'weights'")
    }
  }
  
  if (is(stat, "MCSimulationVector")) {
    if (!identical(attr(stat, "n"), attr(stat, "keyList")[[1]])) {
      stop("attribute n of stat does coincide with attribute keyList of stat")
    }
    
    if (!identical(data$key, attr(stat, "keyList")[[2]])) {
      stop("stat is given for the wrong parametric family or for wrong parameters")
    }
    
    if (!identical(intervalSystem$intervalSystem, attr(stat, "keyList")[[3]])) {
      stop("stat is given for the wrong intervalSystem")
    }
    
    .saveMatrix(stat = stat, options = options, n = data$n)
    
    n <- attr(stat, "n")
    keyList <- attr(stat, "keyList")
    save <- attr(stat, "save")
    if ((data$type == 100L || data$type == 102L) && !all(intervalSystem$lengths %in% attr(stat, "lengths"))) {
      stop("stat does not contain all required lengths")
    }
    stat <- stat[.inOrdered(attr(stat, "lengths"), intervalSystem$lengths), , drop = FALSE]
    attr(stat, "n") <- n
    attr(stat, "keyList") <- keyList
    attr(stat, "save") <- save
  } else {
    stop("stat must be an object of class 'MCSimulationVector' or 'MCSimulationMaximum'")
  }
  
  stat
}

.getVectorStat <- function(stat, data, intervalSystem, penalty = NULL, r = NULL, options = NULL, ...) {
  if (is.null(stat)) {
    stat <- .loadVector(data = data, intervalSystem = intervalSystem, penalty = penalty,
                        r = r, options = options, ...)
  } 

  if (is(stat, "MCSimulationVector")) {
    stat <- .getStat(stat = stat, data = data, intervalSystem = intervalSystem,
                     r = r, options = options, ...)
    
    n <- attr(stat, "n")
    dataKey <- attr(stat, "keyList")[[2]]
    intervalSystemList <- attr(stat, "keyList")[[3]]
    save <- attr(stat, "save")
    
    stat[stat == -Inf] <- 0

    stat <- switch(penalty,
                   "sqrt" = .colMax(sqrt(2 * stat) - 
                                      sqrt(2 * log(exp(1) * n / (intervalSystem$lengths + data$penaltyShift)))),
                   "log" = .colMax(stat - log(exp(1) * n / (intervalSystem$lengths + data$penaltyShift))),
                   "none" = .colMax(stat),
                   stop("unexpected error: unknown penalty"))
    
    class(stat) <- c("MCSimulationMaximum", class(stat))
    attr(stat, "keyList") <- list(n, dataKey, intervalSystemList, intervalSystem$lengths, penalty)
    attr(stat, "key") <- digest::digest(attr(stat, "keyList"))
    attr(stat, "n") <- n
    attr(stat, "lengths") <- intervalSystem$lengths
    attr(stat, "save") <- save
  }
  
  if (is(stat, "MCSimulationMaximum")) {
    if (!identical(attr(stat, "n"), attr(stat, "keyList")[[1]])) {
      stop("attribute n of stat does coincide with attribute keyList of stat")
    }
    
    if (!identical(data$key, attr(stat, "keyList")[[2]])) {
      stop("stat is given for the wrong parametric family or for wrong parameters")
    }
    
    if (!identical(intervalSystem$intervalSystem, attr(stat, "keyList")[[3]])) {
      stop("stat is given for the wrong intervalSystem")
    }
    
    if (!identical(intervalSystem$lengths, attr(stat, "keyList")[[4]])) {
      stop("stat is given for the wrong lengths")
    }
    
    if (!identical(penalty, attr(stat, "keyList")[[5]])) {
      stop("stat is given for the wrong penalty")
    }
    
    .saveVector(stat = stat, options = options, n = data$n)
  } else {
    stop("stat must be an object of class 'MCSimulationVector' or 'MCSimulationMaximum'")
  }

  stat
}

.critValVector <- function(q, data, intervalSystem, penalty) {
  if (length(q) != length(intervalSystem$lengths)) {
    if (!is.null(attr(q, "n"))) {
      n <- attr(q, "n")
      helpData <- data$MC(data, attr(q, "n"))
      attrLengths <- .intervalSystem(intervalSystem$intervalSystem, lengths = NULL, data = helpData)$lengths
    } else {
      n <- data$n
    }
    
    if (length(q) == data$n) {
      q <- q[.inOrdered(1:data$n, intervalSystem$lengths)]
      attr(q, "n") <- n
    } else if (!is.null(attr(q, "n")) && length(q) == attr(q, "n")) {
      q <- q[.inOrdered(1:attr(q, "n"), intervalSystem$lengths)]
      attr(q, "n") <- n
    } else if (!is.null(attr(q, "n")) && length(q) == length(attrLengths)) {
      q <- q[.inOrdered(attrLengths, intervalSystem$lengths)]
      attr(q, "n") <- n
    } else if (length(q) == length(intervalSystem$possibleLengths)) {
      q <- q[.inOrdered(intervalSystem$possibleLengths, intervalSystem$lengths)]
      attr(q, "n") <- n
    } else if (length(q) == 1L) {
      if (penalty == "sqrt") {
        q <- (q + sqrt(2 * log(exp(1) * n / (intervalSystem$lengths + data$penaltyShift))))^2 / 2
      } else if (penalty == "log") {
        q <- q + log(exp(1) * n / (intervalSystem$lengths + data$penaltyShift))
      } else if (penalty == "none") {
        q <- rep(q, length(intervalSystem$lengths))
      } else if (penalty == "weights") {
        stop("q can not be a single value for penalty 'weights'")
      }
      attr(q, "n") <- n
    } else {
      if (is.null(attr(q, "n"))) {
        stop("vector 'q' has wrong length (", length(q), "), it must be either of length 1 or equal ",
             "to the length of lengths (",length(intervalSystem$lengths), ") ", 
             "or equal to the number of observations 'n' (", data$n, ") or ",
             "equal to the length of all possible lengths for the interval system and parametric family (",
             length(intervalSystem$possibleLengths), "). ",
             "Further options are possible if 'q' has an attr 'n' giving the number of observations for which q ",
             "was computed. For more details on all options see also the section ", 
             "'Computation of critical values / global quantile' in the documentation.")
      } else {
        stop("vector 'q' has wrong length (", length(q), "), it must be either of length 1 or equal ",
             "to the length of lengths (",length(intervalSystem$lengths), ") ", 
             "or equal to the number of observations 'n' (", data$n, ") or ",
             "equal to the length of all possible lengths for the interval system and parametric family (",
             length(intervalSystem$possibleLengths), ") or ",
             "equal to the attribute 'n' of 'q' (", attr(q, "n"), ") ",
             "or equal to the length of all possible lengths for the interval system and parametric family ",
             "for the number of observations given in the attribute 'n' of 'q' (", length(attrLengths), "). ",
             "For more details on all options see also the section ", 
             "'Computation of critical values / global quantile' in the documentation.")
      }
    }
  } else {
    if (is.null(attr(q, "n"))) {
      attr(q, "n") <- data$n
    }
  }
  
  q
}

# returns a vector of length nrow(stat) (will be intervalSystem$lengths)
# contaning the scale dependent critical values
.critValWeights <- function(stat, alpha, weights) {
  if (nrow(stat) == 0) {
    return(numeric(0))
  }
  
  if (is.null(weights)) {
    weights <- rep(1 / nrow(stat), nrow(stat))
  } else {
    if (!is.numeric(weights) || length(weights) != nrow(stat) || any(!is.finite(weights)) || any(weights <= 0)) {
      stop("weights must be a finite positive numeric vector of length ", nrow(stat),
           " (equal to length of lengths)")
    }
    
    if (sum(weights) != 1) {
      weights <- weights / sum(weights)
    }
  }
  
  .criticalValuesWeights(stat, weights, alpha)
}

.loadVector <- function(data, intervalSystem, penalty, r, options, ...) {
  if (is.null(options$load$RDSfile)) {
    stat <- NULL
    if (is.null(list(...)$rand.gen)) {
      stat <- .loadVectorSingle(data = data, intervalSystem = intervalSystem, penalty = penalty, r = r,
                                options = options, increased = FALSE, ...)
      
      if (is.null(stat)) {
        stat <- .loadMatrixSingle(data = data, intervalSystem = intervalSystem, r = r, options = options,
                                  increased = FALSE, ...)
      }
      
      if (is.null(stat)) {
        stat <- .loadVectorSingle(data = data, intervalSystem = intervalSystem, penalty = penalty, r = r,
                                  options = options, increased = TRUE, ...)
      }
      
      if (is.null(stat)) {
        stat <- .loadMatrixSingle(data = data, intervalSystem = intervalSystem, r = r, options = options,
                                  increased = TRUE, ...)
      }
    }
    
    if (is.null(stat)) {
      if (options$simulation == "vector") {
        data <- data$MC(data, data$n)
        intervalSystem <- .intervalSystem(intervalSystem = intervalSystem$intervalSystem,
                                          lengths = intervalSystem$lengths, data = data)
        data$signal <- data$generateSignal(data = data, intervalSystem = intervalSystem)
        stat <- .monteCarloSimulation(r = r, output = "value", data = data, intervalSystem = intervalSystem,
                                      penalty = penalty, ...)
      } else if (options$simulation == "vectorIncreased") {
        data <- data$MC(data, data$nq)
        intervalSystem <- .intervalSystem(intervalSystem = intervalSystem$intervalSystem,
                                          lengths = intervalSystem$lengths, data = data)
        data$signal <- data$generateSignal(data = data, intervalSystem = intervalSystem)
        stat <- .monteCarloSimulation(r = r, output = "value", data = data, intervalSystem = intervalSystem,
                                      penalty = penalty, ...)
      } else if (options$simulation == "matrix") {
        data <- data$MC(data, data$n)
        lengths <- intervalSystem$lengths
        intervalSystem <- .intervalSystem(intervalSystem = intervalSystem$intervalSystem,
                                          lengths = NULL, data = data)
        if (!all(lengths %in% intervalSystem$lengths)) {
          data$save <- FALSE
          intervalSystem <- .intervalSystem(intervalSystem = intervalSystem$intervalSystem,
                                            lengths = lengths, data = data)
        }
        data$signal <- data$generateSignal(data = data, intervalSystem = intervalSystem)
        stat <- .monteCarloSimulation(r = r, output = "vector", data = data, intervalSystem = intervalSystem,
                                      penalty = "none", ...)
      } else if (options$simulation == "matrixIncreased") {
        data <- data$MC(data, data$nq)
        lengths <- intervalSystem$lengths
        intervalSystem <- .intervalSystem(intervalSystem = intervalSystem$intervalSystem,
                                          lengths = NULL, data = data)
        if (!all(lengths %in% intervalSystem$lengths)) {
          data$save <- FALSE
          intervalSystem <- .intervalSystem(intervalSystem = intervalSystem$intervalSystem,
                                            lengths = lengths, data = data)
        }
        data$signal <- data$generateSignal(data = data, intervalSystem = intervalSystem)
        stat <- .monteCarloSimulation(r = r, output = "vector", data = data, intervalSystem = intervalSystem,
                                      penalty = "none", ...)
      }
    }
  } else {
    stat <- readRDS(options$load$RDSfile)
  }
  
  stat
}

.loadMatrix <- function(data, intervalSystem, r, options, increased = FALSE, ...) {
  if (is.null(options$load$RDSfile)) {
    stat <- NULL
    if (is.null(list(...)$rand.gen)) {
      stat <- .loadMatrixSingle(data = data, intervalSystem = intervalSystem, r = r, options = options,
                                increased = FALSE, ...)
      
      if (is.null(stat)) {
        stat <- .loadMatrixSingle(data = data, intervalSystem = intervalSystem, r = r, options = options,
                                  increased = TRUE, ...)
      }
    }
    
    if (is.null(stat)) {
      if (options$simulation == "matrix") {
        data <- data$MC(data, data$n)
        lengths <- intervalSystem$lengths
        intervalSystem <- .intervalSystem(intervalSystem = intervalSystem$intervalSystem,
                                          lengths = NULL, data = data)
        if (!all(lengths %in% intervalSystem$lengths)) {
          data$save <- FALSE
          intervalSystem <- .intervalSystem(intervalSystem = intervalSystem$intervalSystem,
                                            lengths = lengths, data = data)
        }
        data$signal <- data$generateSignal(data = data, intervalSystem = intervalSystem)
        stat <- .monteCarloSimulation(r = r, output = "vector", data = data, intervalSystem = intervalSystem,
                                      penalty = "none", ...)
      } else if (options$simulation == "matrixIncreased") {
        data <- data$MC(data, data$nq)
        lengths <- intervalSystem$lengths
        intervalSystem <- .intervalSystem(intervalSystem = intervalSystem$intervalSystem,
                                          lengths = NULL, data = data)
        if (!all(lengths %in% intervalSystem$lengths)) {
          data$save <- FALSE
          intervalSystem <- .intervalSystem(intervalSystem = intervalSystem$intervalSystem,
                                            lengths = lengths, data = data)
        }
        data$signal <- data$generateSignal(data = data, intervalSystem = intervalSystem)
        stat <- .monteCarloSimulation(r = r, output = "vector", data = data, intervalSystem = intervalSystem,
                                      penalty = "none", ...)
      }
    }
  } else {
    stat <- readRDS(options$load$RDSfile)
  }
  
  stat
}

.loadMatrixSingle <- function(data, intervalSystem, r, options, increased = FALSE, ...) {
  stat <- NULL
  
  if (increased) {
    keyList <- list(data$nq, data$key, intervalSystem$intervalSystem)
    key <- digest::digest(keyList)
    look <- "matrixIncreased"
  } else {
    keyList <- list(data$n, data$key, intervalSystem$intervalSystem)
    key <- digest::digest(keyList)
    look <- "matrix"
  }
  
  if (!is.null(options$load)) {
    if (!is.null(options$load$workspace)) {
      if (look %in% options$load$workspace) {
        stat <- .loadWorkspace(key = key, r = r, envir = options$envir)
      }
    }
    
    if (is.null(stat) && !is.null(options$load$fileSystem)) {
      if (look %in% options$load$fileSystem) {
        stat <- .loadRcache(key = keyList, r = r, dirs = options$dirs)
      }
    }
    
    if (is.null(stat) && increased && !is.null(options$load$package) && options$load$package) {
      stat <- .loadPackage(data, intervalSystem, r)
    }
  }
  
  stat
}

.loadVectorSingle <- function(data, intervalSystem, penalty, r, options, increased = FALSE, ...) {
  stat <- NULL
  
  if (increased) {
    keyList <- list(data$nq, data$key, intervalSystem$intervalSystem, intervalSystem$lengths, penalty)
    key <- digest::digest(keyList)
    look <- "vectorIncreased"
  } else {
    keyList <- list(data$n, data$key, intervalSystem$intervalSystem, intervalSystem$lengths, penalty)
    key <- digest::digest(keyList)
    look <- "vector"
  }
  
  if (!is.null(options$load)) {
    if (!is.null(options$load$workspace)) {
      if (look %in% options$load$workspace) {
        stat <- .loadWorkspace(key = key, r = r, envir = options$envir)
      }
    }
    
    if (is.null(stat) && !is.null(options$load$fileSystem)) {
      if (look %in% options$load$fileSystem) {
        stat <- .loadRcache(key = keyList, r = r, dirs = options$dirs)
      }
    }
  }
  
  stat
}

.saveMatrix <- function(stat, options, n) {
  if (!is.null(options$save)) {
    if (!is.null(options$save$workspace) && attr(stat, "save")) {
      if (attr(stat, "n") == n) {
        if ("matrix" %in% options$save$workspace) {
          .saveWorkspace(key = attr(stat, "key"), stat = stat, r = ncol(stat), envir = options$envir)
        }
      } else {
        if ("matrixIncreased" %in% options$save$workspace) {
          .saveWorkspace(key = attr(stat, "key"), stat = stat, r = ncol(stat), envir = options$envir)
        }
      }
    }
    
    if (!is.null(options$save$fileSystem) && attr(stat, "save")) {
      if (attr(stat, "n") == n) {
        if ("matrix" %in% options$save$fileSystem) {
          .saveRcache(key = attr(stat, "keyList"), stat = stat, r = ncol(stat), dirs = options$dirs)
        }
      } else {
        if ("matrixIncreased" %in% options$save$fileSystem) {
          .saveRcache(key = attr(stat, "keyList"), stat = stat, r = ncol(stat), dirs = options$dirs)
        }
      }
    }
    
    if (!is.null(options$save$RDSfile) && options$save$RDSfile[2] != "") {
      saveRDS(stat, file = options$save$RDSfile[2])
    }
    
    if (!is.null(options$save$variable) && options$save$variable[2] != "") {
      assign(options$save$variable[2], stat, pos = options$envir)
    }
  }
}

.saveVector <- function(stat, options, n) {
  if (!is.null(options$save)) {
    if (!is.null(options$save$workspace) && attr(stat, "save")) {
      if (attr(stat, "n") == n) {
        if ("vector" %in% options$save$workspace) {
          .saveWorkspace(key = attr(stat, "key"), stat = stat, r = length(stat), envir = options$envir)
        }
      } else {
        if ("vectorIncreased" %in% options$save$workspace) {
          .saveWorkspace(key = attr(stat, "key"), stat = stat, r = length(stat), envir = options$envir)
        }
      }
    }
    
    if (!is.null(options$save$fileSystem) && attr(stat, "save")) {
      if (attr(stat, "n") == n) {
        if ("vector" %in% options$save$fileSystem) {
          .saveRcache(key = attr(stat, "keyList"), stat = stat, r = length(stat), dirs = options$dirs)
        }
      } else {
        if ("vectorIncreased" %in% options$save$fileSystem) {
          .saveRcache(key = attr(stat, "keyList"), stat = stat, r = length(stat), dirs = options$dirs)
        }
      }
    }
    
    if (!is.null(options$save$RDSfile) && options$save$RDSfile[1] != "") {
      saveRDS(stat, file = options$save$RDSfile[1])
    }
    
    if (!is.null(options$save$variable) && options$save$variable[1] != "") {
      assign(options$save$variable[1], stat, pos = options$envir)
    }
  }
}

.loadWorkspace <- function(key, r, envir) {
  stat <- NULL
  if (exists("critValStepRTab", where = envir, inherits = FALSE)) {
    tab <- get("critValStepRTab", pos = envir, inherits = FALSE)
    
    index <- which(tab$key == key)
    if (length(index) == 1 && tab$r[index] >= r) {
      stat <- tab$stat[[index]]
    }
  }
  stat
}

.saveWorkspace <- function(key, stat, r, envir) {
  if (exists("critValStepRTab", where = envir, inherits = FALSE)) {
    tab <- get("critValStepRTab", pos = envir, inherits = FALSE)
    
    if (all(key != tab$key)) {
      index <- length(tab$key) + 1L
      tab$key[index] <- key
      tab$r[index] <- r
      tab$stat[[index]] <- stat
    } else {
      index <- which(key == tab$key)
      if (r > tab$r[index]) {
        tab$key[index] <- key
        tab$r[index] <- r
        tab$stat[[index]] <- stat
      }
    }
  } else {
    tab <- list(key = key, r = r, stat = list(stat))
  }
  assign("critValStepRTab", tab, pos = envir, inherits = FALSE)
}

.loadRcache <- function(key, r, dirs) {
  stat <- R.cache::loadCache(key = key, dirs = dirs)
  if (!is.null(stat)) {
    if (is.matrix(stat)) {
      if (ncol(stat) < r) {
        stat <- NULL
      }
    } else {
      if (length(stat) < r) {
        stat <- NULL
      }
    }
  }
  stat
}

.saveRcache <- function(key, stat, r, dirs) {
  old <- R.cache::loadCache(key = key, dirs = dirs)
  if (is.null(old)) {
    R.cache::saveCache(stat, key = key, dirs = dirs)
  } else {
    if (is.matrix(old)) {
      if (ncol(old) < r) {
        R.cache::saveCache(stat, key = key, dirs = dirs)
      }
    } else {
      if (length(old) < r) {
        R.cache::saveCache(stat, key = key, dirs = dirs)
      }
    }
  }
}

.loadPackage <- function(data, intervalSystem, r) {
  stat <- NULL
  if (r <= 1e4) {
    path <- system.file("extdata", sep = "", package = "stepRdata")
    path <- file.path(path, data$family, paste(intervalSystem$intervalSystem, data$nq, ".rds", sep = ""))
    
    if (file.exists(path)) {
      stat <- readRDS(path)
      
      # to be safe if results of digest change
      if (data$family == "gauss") {
        attr(stat, "keyList")[[2]] <- digest::sha1(list("gauss"), digits = 6)
      }
      
      if (data$family == "hsmuce") {
        attr(stat, "keyList")[[2]] <- digest::sha1(list("hsmuce"), digits = 6)
      }
      attr(stat, "key") <- digest::digest(attr(stat, "keyList"))
    }
  }
  
  stat
}
