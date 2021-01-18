.intervalSystem <- function(intervalSystem, lengths = NULL, data) {
  if (is.null(intervalSystem)) {
    intervalSystem <- data$defaultIntervalSystem
  }
  intervalSystem <- match.arg(intervalSystem, c("all", "dyaLen", "dyaPar"))


  ret <- list(intervalSystem = intervalSystem)

  possibleLengthsIntervalSystem <- switch(intervalSystem,
                                          all = 1:data$n,
                                          dyaLen = 2^(0:as.integer(floor(log2(data$n)) + 1e-6)),
                                          dyaPar = 2^(0:as.integer(floor(log2(data$n)) + 1e-6))
  )
  ret$possibleLengths <- data$possibleLengths[.inOrdered(data$possibleLengths, possibleLengthsIntervalSystem)]

  if (is.null(lengths)) {
    lengths <- data$defaultLengths[.inOrdered(data$defaultLengths, ret$possibleLengths)]
  } else {
    if (!is.numeric(lengths) || any(!is.finite(lengths))) {
      stop("lengths must be an integer vector containing finite values")
    }
    
    if (any(!is.integer(lengths))) {
      lengths <- as.integer(lengths + 1e-6)
    }
    
    if (is.unsorted(lengths, strictly = TRUE)) {
      lengths <- sort(lengths)
      if (is.unsorted(lengths, strictly = TRUE)) {
        warning("lengths contains duplicated values, they will be removed")
        lengths <- unique(lengths)
      }
    }

    if (any(!(.inOrdered(lengths, ret$possibleLengths)))) {
      wrongIntervalSystem <- !(.inOrdered(lengths, possibleLengthsIntervalSystem))
      wrongData <- !(.inOrdered(lengths, data$possibleLengths))
      if (any(wrongIntervalSystem) && any(wrongData)) {
        stop("argument 'lengths' contains inappropriate values.",
             " The following lengths are not possible for ", data$n, " observations and ",
             "interval system '", intervalSystem, "': ", paste(lengths[wrongIntervalSystem], collapse = ', '),
             ". ", "And the following lengths are not possible for parametric family '", data$family, "': ",
             paste(lengths[wrongData], collapse = ', '), ". ", 
             "Please see also the documentation for possible choices.")
      } else if (any(wrongIntervalSystem)) {
        stop("argument 'lengths' contains inappropriate values.",
             " The following lengths are not possible for ", data$n, " observations and ",
             "interval system '", intervalSystem, "': ", paste(lengths[wrongIntervalSystem], collapse = ', '),
             ". ", "Please see also the documentation for possible choices.")
      } else if (any(wrongData)) {
        stop("argument 'lengths' contains inappropriate values.",
             " The following lengths are not possible for ", data$n, " observations and ",
             "parametric family '", data$family, "': ", paste(lengths[wrongData], collapse = ', '), ". ",
             "Please see also the documentation for possible choices.")
      }
    }
  }

  ret$lengths <- lengths
  
  if (intervalSystem == "all") {
    lengths <- logical(data$n)
    lengths[ret$lengths] <- TRUE
    if (all(lengths)) {
      ret$type = 0L
    } else {
      ret$type = 1L
      ret$argumentsList = list(lengths = lengths)
    }
  } else if (intervalSystem == "dyaLen") {
    if (length(lengths) == as.integer(floor(log2(data$n)) + 1e-6) + 1) {
      ret$type = 10L
    } else {
      ret$type = 11L
      ret$argumentsList = list(lengths = as.integer(ret$lengths))
    }    
  } else if (intervalSystem == "dyaPar") {
    if (length(lengths) == as.integer(floor(log2(data$n)) + 1e-6) + 1) {
      ret$type = 20L
    } else {
      ret$type = 21L
      ret$argumentsList = list(lengths = as.integer(ret$lengths))
    }    
  }
  ret
}
