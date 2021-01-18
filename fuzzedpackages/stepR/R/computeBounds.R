computeBounds <- function(y, q = NULL, alpha = NULL, family = NULL, intervalSystem = NULL, lengths = NULL, ...) {
  if (missing(y)) {
    stop("argument 'y' must be given")
  }
  
  .RemoveAdditionalArgsPF <- function(family, y, n, ..., penalty, alpha, stat, r, weights, options,
                                      seed, rand.gen, messages)
    .parametricFamily(family = family, y = y, n = n, ...)
  data <- .RemoveAdditionalArgsPF(family = family, y = y, n = length(y), ...)
  intervalSystem <- .intervalSystem(intervalSystem = intervalSystem, lengths = lengths, data = data)
  
  .RemoveAdditionalArgsCV <- function(q, alpha, data, output, intervalSystem, ..., nq,
                                      sd, covariances, correlations, filter)
    .critVal(q = q, alpha = alpha, data = data, output = output, intervalSystem = intervalSystem, ...)
  q <- .RemoveAdditionalArgsCV(q = q, alpha = alpha, data = data, output = "vector",
                               intervalSystem = intervalSystem, ...)
  
  criticalValues <- rep(Inf, data$n)
  criticalValues[intervalSystem$lengths] <- q
  
  as.data.frame(.callRoutines(observations = data$y, routineType = 2L,
                              argumentsListRoutine = list(q = criticalValues), 
                              dataType = data$type, argumentsListData = data$argumentsList,
                              intervalSystemType = intervalSystem$type,
                              argumentsListIntervalSystem = intervalSystem$argumentsList))
}
