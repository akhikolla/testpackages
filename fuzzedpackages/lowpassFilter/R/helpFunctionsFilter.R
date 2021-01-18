getSignalJump <- function(t, cp, leftValue, rightValue) {
  ifelse(t < cp, rep(leftValue, length(t)), rep(rightValue, length(t)))
}

getConvolutionJump <- function(t, cp, leftValue, rightValue, filter, truncated = TRUE) {
  if (truncated) {
    ret <- leftValue * (1 - filter$truncatedStepfun(t - cp)) + rightValue * filter$truncatedStepfun(t - cp)
  } else {
    ret <- leftValue * (1 - filter$stepfun(t - cp)) + rightValue * filter$stepfun(t - cp)
  }
  ret
}

getSignalPeak <- function(t, cp1, cp2, value, leftValue, rightValue) {
  ifelse(t < cp1, rep(leftValue, length(t)), ifelse(t < cp2, rep(value, length(t)), rep(rightValue, length(t))))
}

getConvolutionPeak <- function(t, cp1, cp2, value, leftValue, rightValue, filter, truncated = TRUE) {
  if (truncated) {
    ret <- leftValue * (1 - filter$truncatedStepfun(t - cp1)) + 
      value * (filter$truncatedStepfun(t - cp1) - filter$truncatedStepfun(t - cp2)) + 
      rightValue * filter$truncatedStepfun(t - cp2)
  } else {
    ret <- leftValue * (1 - filter$stepfun(t - cp1)) + 
      value * (filter$stepfun(t - cp1) - filter$stepfun(t - cp2)) + 
      rightValue * filter$stepfun(t - cp2)
  }
  ret
}

getConvolution <- function(t, stepfun, filter, truncated = TRUE) {
  sapply(t, function(s) {
    arguments <- s - c(stepfun$leftEnd[1], stepfun$rightEnd)
    if (truncated) {
      steps <- filter$truncatedStepfun(arguments)
    } else {
      steps <- filter$stepfun(arguments)
    }
    sum(stepfun$value * (c(1, steps[-c(1, length(steps))]) - c(steps[-c(1, length(steps))], 0)))
  })
}
