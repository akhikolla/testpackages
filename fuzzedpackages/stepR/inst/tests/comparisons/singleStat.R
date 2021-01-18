###
# single statistic, has to be implemented for all parametric families
# singleStat <- function(y, value, li, ri, ...)
# return single numeric
###

# family gauss: independent homogeneous gaussian observations
singleStatGauss <- function(y, value, li, ri, sd, ...) {
  length(y[li:ri]) * (mean(y[li:ri]) - value)^2 / 2 / sd^2
}

# family mDependentPS: partial sum statistic for dependent homogeneous gaussian observations
singleStatmDependentPS <- function(y, value, li, ri, covariances, ...) {
  len <- ri - li + 1
  m <- min(len, length(covariances) - 1L)
  varianceSum <- numeric(1)
  if (m == 0) {
    varianceSum <- len * covariances[1]
  } else {
    varianceSum <- len * covariances[1] + 2 * sum((len - 1:m) * covariances[2:(m + 1)])
  }
  sum(y[li:ri] - value)^2 / varianceSum / 2
}

# family jsmurf: m-dependent homogeneous gaussian observations using jsmurf idea
singleStatJsmurf <- function(y, value, li, ri, sd, filter, ...) {
  if (ri - li < filter$len) {
    return(-Inf)
  }
  li <- li + filter$len
  length(y[li:ri]) * (mean(y[li:ri]) - value)^2 / 2 / sd^2
}

# family jsmurftPS: partial sum statistic for dependent homogeneous gaussian observations using jsmurf idea
singleStatJsmurfPS <- function(y, value, li, ri, sd, filter, ...) {
  if (ri - li < filter$len) {
    return(-Inf)
  }
  li <- li + filter$len
  len <- ri - li + 1
  covariances <- filter$acf * sd^2
  m <- min(len, length(covariances) - 1L)
  varianceSum <- numeric(1)
  if (m == 0) {
    varianceSum <- len * covariances[1]
  } else {
    varianceSum <- len * covariances[1] + 2 * sum((len - 1:m) * covariances[2:(m + 1)])
  }
  (sum(y[li:ri]) - value * len)^2 / varianceSum / 2
}

# family jsmurfLR: likelihood ratio test statistic for dependent homogeneous gaussian observations
#                  using jsmurf idea
singleStatJsmurfLR <- function(y, value, li, ri, sd, filter, ...) {
  if (ri - li < filter$len) {
    return(-Inf)
  }
  li <- li + filter$len
  len <- ri - li + 1
  covariances <- filter$acf * sd^2
  m <- min(len, length(covariances) - 1L)
  
  if (len == 1) {
    return((y[li] - value)^2 / covariances[1] / 2)
  }
  
  A <- matrix(0, len, len)
  for (i in 1:(len - 1)) {
    A[i, i] <- covariances[1]
    A[i, i + 1:min(m, len - i)] <- covariances[2:min(m + 1, len - i + 1)]
    A[i + 1:min(m, len - i), i] <- covariances[2:min(m + 1, len - i + 1)]
  }
  A[len, len] <- covariances[1]  
  
  inverseAone <- solve(A, rep(1, len))
  sum((y[li:ri] - value) * inverseAone)^2 / sum(inverseAone) / 2
}

# family hsmuce: independent heterogeneous gaussian observations
singleStatHsmuce <- function(y, value, li, ri, ...) {
  length(y[li:ri]) * (mean(y[li:ri]) - value)^2 / stats::var(y[li:ri]) / 2
}

# family hjsmurf: m-dependent heterogeneous gaussian observations using jsmurf idea
singleStatHjsmurf <- function(y, value, li, ri, filter, ...) {
  if (ri - li <= filter$len) {
    return(-Inf)
  }
  li <- li + filter$len
  length(y[li:ri]) * (mean(y[li:ri]) - value)^2 / stats::var(y[li:ri]) / 2
}

# family hjsmurfSPS: partial sum statistic for dependent heterogeneous gaussian observations using jsmurf idea
singleStatHjsmurfSPS <- function(y, value, li, ri, filter, ...) {
  if (ri - li <= filter$len) {
    return(-Inf)
  }
  li <- li + filter$len
  len <- ri - li + 1
  correlations <- filter$acf[-1]
  m <- min(len - 1, length(correlations))
  varianceSum <- numeric(1)
  if (m == 0) {
    varianceSum <- len
  } else {
    varianceSum <- len + 2 * sum((len - 1:m) * correlations[1:m])
  }
  estVar <- stats::var(y[li:ri]) * (len - 1) / len / (1 - varianceSum / len^2)
  
  (sum(y[li:ri]) - value * len)^2 / varianceSum / 2 / estVar
}

# family hsjmurfLR: likelihood ratio test statistic for dependent homogeneous gaussian observations
#                   using jsmurf idea
singleStatHjsmurfLR <- function(y, value, li, ri, filter, ...) {
  if (ri - li <= filter$len) {
    return(-Inf)
  }
  li <- li + filter$len
  len <- ri - li + 1
  correlations <- filter$acf[-1]
  m <- min(len, length(correlations))
  
  if (len == 1) {
    return(-Inf)
  }
  
  A <- matrix(0, len, len)
  for (i in 1:(len - 1)) {
    A[i, i] <- 1
    A[i, i + 1:min(m, len - i)] <- correlations[1:min(m, len - i)]
    A[i + 1:min(m, len - i), i] <- correlations[1:min(m, len - i)]
  }
  A[len, len] <- 1  
  
  inverseAyMean <- solve(A, y[li:ri] - mean(y[li:ri]))
  inverseAyValue <- solve(A, y[li:ri] - value)
  sum((y[li:ri] - value) * inverseAyValue) / sum((y[li:ri] - mean(y[li:ri])) * inverseAyMean) / 2
}

#############
# local tests

singleStat2Param <- function(obs, time, filter, left, right, leftValue, rightValue,
                                    leftVar, rightVar, cp, ...) {
  if (is.na(leftVar) || is.na(rightVar)) {
    return(NA)
  }
  
  Fleft <- filter$truncatedStepfun(time - left)
  Fright <- filter$truncatedStepfun(time - right)
  obs2 <- obs
  obs <- (obs - leftValue * (1 - Fleft) - rightValue * Fright)
  v <- Fleft - Fright
  sumv2 <- sum(v^2)
  est <- sum(obs * v) / sumv2
  
  Fleft <- outer(time, time, function(i, j) filter$acAntiderivative(pmin(i, j) - left, abs(j - i)))
  Fright <- outer(time, time, function(i, j) filter$acAntiderivative(pmin(i, j) - right, abs(j - i)))
  cor <- outer(time, time, function(i, j) filter$acfun(abs(j - i)))
  w <- Fleft - Fright
  sigmaLR <- leftVar * (cor - Fleft) + rightVar * Fright
  vv <- outer(seq(along = time), seq(along = time), function(i, j) v[i] * v[j] / sum(v^2))
  diagW <- diag(w)
  matrixDiagW <- matrix(rep(diagW, length(diagW)), length(diagW))
  A <- sum(diag(sigmaLR) * diagW) - sum(vv * sigmaLR * matrixDiagW)
  B <- sum(diag(w) * diagW) - sum(vv * w * matrixDiagW)
  sigmaest2 <- (sum(diagW * (obs - v * est)^2) - A) / B
  
  if (is.na(sigmaest2)) {
    return(NA)
  }
  
  if (sigmaest2 < 0) {
    sigmaest2 <- 0
  }
  
  w <- diag(w)
  sigmaLR <- diag(sigmaLR)
  Tvar <- w * sigmaest2 + sigmaLR
  Test <- (obs - v * est)^2
  
  Fcp <- filter$truncatedStepfun(time - cp)
  obs2 <- (obs2 - leftValue * (1 - Fcp) - rightValue * Fcp)
  Tobs <- obs2^2
  
  Fcp <-  filter$acAntiderivative(time - cp, 0)
  T0 <- leftVar * (1 - Fcp) + rightVar * Fcp
  
  sum(log(T0 / Tvar) + Tobs / T0 - Test / Tvar)
}

singleStatLR <- function(obs, time, filter, left, right, leftValue, rightValue, cp, sd, regu, ...) {
  covariances <- filter$acf * sd^2
  covariances[1] <- covariances[1] + regu * sd^2
  
  len <- length(obs)
  m <- min(len, length(covariances) - 1L)
  
  if (len == 1) {
    A <- matrix(covariances[1], 1, 1)
  } else {
    A <- matrix(0, len, len)
    for (i in 1:(len - 1)) {
      A[i, i] <- covariances[1]
      A[i, i + 1:min(m, len - i)] <- covariances[2:min(m + 1, len - i + 1)]
      A[i + 1:min(m, len - i), i] <- covariances[2:min(m + 1, len - i + 1)]
    }
    A[len, len] <- covariances[1]  
  }
  
  Fleft <- filter$truncatedStepfun(time - left)
  Fright <- filter$truncatedStepfun(time - right)
  v <- Fleft - Fright
  sol <- solve(A, v)
  
  est <- sum((obs - leftValue * (1 - Fleft) - rightValue * Fright) * sol) / sum(v * sol)
  
  convolvedSignal <- lowpassFilter::getConvolutionPeak(time, left, right, est, leftValue, rightValue, filter)
  convolvedNull <- lowpassFilter::getConvolutionJump(time, cp, leftValue, rightValue, filter)
  
  sum((obs - convolvedNull) * solve(A, obs - convolvedNull)) - 
    sum((obs - convolvedSignal) * solve(A, obs - convolvedSignal))
}
