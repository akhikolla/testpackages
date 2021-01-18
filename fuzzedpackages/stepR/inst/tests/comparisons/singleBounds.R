###
# single bounds, has to be implemented for all parametric families
# singleBounds <- function(y, criticalValue, li, ri, ...)
# return single numeric
###

# family gauss: independent homogeneous gaussian observations
singleBoundsGauss <- function(y, criticalValue, li, ri, sd) {
  len <- ri - li + 1
  m <- mean(y[li:ri])
  s <- sd * sqrt(2 * criticalValue / len)
  c(m - s, m + s)
}

# family gaussDependent: partial sum statistic for dependent homogeneous gaussian observations
singleBoundsmDependentPS <- function(y, criticalValue, li, ri, covariances) {
  len <- ri - li + 1
  m <- min(len, length(covariances) - 1L)
  varianceSum <- numeric(1)
  if (m == 0) {
    varianceSum <- len * covariances[1]
  } else {
    varianceSum <- len * covariances[1] + 2 * sum((len - 1:m) * covariances[2:(m + 1)])
  }
  m <- mean(y[li:ri])
  s <- sqrt(2 * varianceSum * criticalValue) / len
  c(m - s, m + s)
}

# family jsmurf: m-dependent homogeneous gaussian observations using jsmurf idea
singleBoundsJsmurf <- function(y, criticalValue, li, ri, sd, filter) {
  if (ri - li < filter$len) {
    return(c(-Inf, Inf))
  }
  li <- li + filter$len
  
  len <- ri - li + 1
  m <- mean(y[li:ri])
  s <- sd * sqrt(2 * criticalValue / len)
  c(m - s, m + s)
}

# family jsmurfPS: partial sum statistic for dependent homogeneous gaussian observations using jsmurf idea
singleBoundsJsmurfPS <- function(y, criticalValue, li, ri, sd, filter) {
  if (ri - li < filter$len) {
    return(c(-Inf, Inf))
  }
  li <- li + filter$len
  
  len <- ri - li + 1
  covariances <- sd^2 * filter$acf
  m <- min(len, length(covariances) - 1L)
  varianceSum <- numeric(1)
  if (m == 0) {
    varianceSum <- len * covariances[1]
  } else {
    varianceSum <- len * covariances[1] + 2 * sum((len - 1:m) * covariances[2:(m + 1)])
  }
  m <- mean(y[li:ri])
  s <- sqrt(criticalValue * 2 * varianceSum) / len
  c(m - s, m + s)
}

# family jsmurfLR: likelihood ratio test statistic for dependent homogeneous gaussian observations
#                  using jsmurf idea
singleBoundsJsmurfLR <- function(y, criticalValue, li, ri, sd, filter) {
  if (ri - li < filter$len) {
    return(c(-Inf, Inf))
  }
  li <- li + filter$len
  
  len <- ri - li + 1
  covariances <- sd^2 * filter$acf
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
  
  inverseAone <- solve(A, rep(1, len))
  m <- sum(y[li:ri] * inverseAone)
  s <- sqrt(criticalValue * 2 * sum(inverseAone))
  c(m - s, m + s) / sum(inverseAone)
}

# family hsmuce: independent heterogeneous gaussian observations
singleBoundsHsmuce <- function(y, criticalValue, li, ri) {
  len <- ri - li + 1
  m <- mean(y[li:ri])
  s <- stats::sd(y[li:ri]) * sqrt(2 * criticalValue / len)
  c(m - s, m + s)
}

# family hjsmurf: m-dependent heterogeneous gaussian observations using jsmurf idea
singleBoundsHjsmurf <- function(y, criticalValue, li, ri, filter) {
  if (ri - li <= filter$len) {
    return(c(-Inf, Inf))
  }
  li <- li + filter$len
  
  len <- ri - li + 1
  m <- mean(y[li:ri])
  s <- stats::sd(y[li:ri]) * sqrt(2 * criticalValue / len)
  c(m - s, m + s)
}

# family jsmurfSPS: studentized partial sum statistic for dependent heterogeneous gaussian observations
#                   using jsmurf idea
singleBoundsHjsmurfSPS <- function(y, criticalValue, li, ri, filter) {
  if (ri - li <= filter$len) {
    return(c(-Inf, Inf))
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
  m <- mean(y[li:ri])
  s <- sqrt(2 * criticalValue * varianceSum * estVar) / len
  c(m - s, m + s)
}

# family hjsmurfLR: likelihood ratio test statistic for dependent heterogeneous gaussian observations
#                   using jsmurf idea
singleBoundsHjsmurfLR <- function(y, criticalValue, li, ri, filter) {
  if (ri - li <= filter$len) {
    return(c(-Inf, Inf))
  }
  li <- li + filter$len
  
  len <- ri - li + 1
  correlations <- filter$acf[-1]
  m <- min(len, length(correlations))
  
  A <- matrix(0, len, len)
  for (i in 1:(len - 1)) {
    A[i, i] <- 1
    A[i, i + 1:min(m, len - i)] <- correlations[1:min(m, len - i)]
    A[i + 1:min(m, len - i), i] <- correlations[1:min(m, len - i)]
  }
  A[len, len] <- 1  
  
  inverseA <- solve(A)
  m <- t(y[li:ri]) %*% inverseA %*% rep(1, len)
  s <- sqrt((t(y[li:ri]) %*% inverseA %*% rep(1, len))^2 - (t(rep(1, len)) %*% inverseA %*% rep(1, len)) * 
              (t(y[li:ri]) %*% inverseA %*% y[li:ri] - 2 * criticalValue * 
                 t(y[li:ri] - mean(y[li:ri])) %*% inverseA %*% (y[li:ri] - mean(y[li:ri]))))
  c(m - s, m + s) / as.numeric(t(rep(1, len)) %*% inverseA %*% rep(1, len))
}
