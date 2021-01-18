.testSmallScales <- function(data, family, lengths = NULL, q, alpha, ...) {
  .RemoveAdditionalArgsPF <- function(family, y, n, ..., penalty, alpha, stat, r, weights, options,
                                      seed, rand.gen, messages)
    .parametricFamily(family = family, y = y, n = n, ...)
  data <- .RemoveAdditionalArgsPF(family = family, y = data, n = length(data), ...)
  
  intervalSystem <- .intervalSystem(intervalSystem = NULL, lengths = lengths, data = data)
  signal <- data$generateSignal(data = data, intervalSystem = intervalSystem)
  
  # input check of q
  .RemoveAdditionalArgsCV <- function(q, alpha, data, output, intervalSystem, ..., nq,
                                      filter, fit, startTime, thresholdLongSegment, localValue,
                                      localVar, regularization, correlations, suppressWarningNoDeconvolution, localList)
    .critVal(q = q, alpha = alpha, data = data, output = output, intervalSystem = intervalSystem, ...)
  signal$critVal <- .RemoveAdditionalArgsCV(q = q, alpha = alpha, data = data, output = "vector",
                                            intervalSystem = intervalSystem, ...)
  
  ret <- .callRoutines(observations = data$y, routineType = 10L, argumentsListRoutine = signal, 
                       dataType = data$type, argumentsListData = data$argumentsList,
                       intervalSystemType = intervalSystem$type,
                       argumentsListIntervalSystem = intervalSystem$argumentsList)

  add <- data.frame(left = c(data$addLeft - data$filter$len, as.integer(ret$left - data$filter$len / 2 + data$tolerance)),
                    right = c(data$addRight, as.integer(ret$right + data$filter$len / 2 - data$tolerance) + 1L),
                    noDeconvolution = c(data$noDeconvolution, ret$noDeconvolution))
  add$left[add$left < 1L] <- 1L
  add <- add[order(add$left), ]
  addLeft <- add$left
  addRight <- add$right
  noDeconvolution <- add$noDeconvolution
  
  i <- 2
  while (i <= length(addLeft)) {
    if (addLeft[i] <= addRight[i - 1]) {
      if (i + 1L <= length(noDeconvolution)) {
        if (1L <= i - 2L) {
          noDeconvolution <- c(noDeconvolution[1:(i - 2)], TRUE,
                               noDeconvolution[(i + 1):length(noDeconvolution)])
        } else {
          noDeconvolution <- c(TRUE, noDeconvolution[(i + 1):length(noDeconvolution)])
        }
      } else {
        if (1L <= i - 2L) {
          noDeconvolution <- c(noDeconvolution[1:(i - 2)], TRUE)
        } else {
          noDeconvolution <- TRUE
        }
      }
      if (i + 1L <= length(addLeft)) {
        addLeft <- c(addLeft[1:(i - 1)], addLeft[(i + 1):length(addLeft)])
      } else {
        addLeft <- addLeft[1:(i - 1)]
      }
      if (1L <= i - 2L) {
        addRight <- c(addRight[1:(i - 2)], addRight[i:length(addRight)])
      } else {
        addRight <- addRight[i:length(addRight)]
      }
    } else {
      i <- i + 1
    }
  }
  
  jumps <- as.integer((data$jumps - data$startTime) * data$filter$sr + 1e-6)
  
  included <- logical(length(jumps))
  jumpsD <- as.integer((data$jumpsD - data$startTime) * data$filter$sr + 1e-6)
  for (i in seq(along = addLeft)) {
    included[which(jumpsD >= addLeft[i] - data$filter$len & jumpsD <= addRight[i] + data$filter$len)] <- TRUE
  }
  jumps <- jumps[!included]
  
  list(jumps = jumps, addLeft = addLeft, addRight = addRight, noDeconvolution = noDeconvolution, data = data,
       q = signal$critVal)
}
