# Number of lags and directions
#  1. Directions with all lags
#  2. Lags with all directions
#  3. Combinations of 1 and 2

lagdir_rect <- function(L, thres = 0.8) {
  value <- 0
  lag <- 0
  dir <- 0
  sumsc <- colSums(L)
  sumsr <- rowSums(L)
  while (value < thres) {
    dir <- dir + 1
    value <- sum(sumsc[1:dir])/sum(sumsc)
  }
  alllag <- dir
  value <- 0
  while (value < thres) {
    lag <- lag + 1    
    value <- sum(sumsr[1:lag])/sum(sumsr)
  }
  alldir <- lag
  value <- sum(L[1:lag, 1:dir])
  # If this value > threshold, then result is achieved already
  if (value < thres) {
    ks <- Inf
    lag.old <- lag
    dir.old <- dir
    for (i in dir.old:ncol(L)) { #If 2*n works and 1*n not, 2*(n-a), a=1, 2... might still work!
      if ((i * lag.old) <= ks) {
        dir.new <- i
        value1 <- sum(L[1:lag.old, 1:dir.new])
        lag.new <- lag.old
        while ((lag.new < nrow(L)) && (value1 < thres)) {
          lag.new <- lag.new + 1
          value1 <- sum(L[1:lag.new, 1:dir.new])
        }
        ks.new <- dir.new*lag.new #k*s rectangle
        if (ks.new < ks) {
          value2 <- value1
          ks <- ks.new
          dir <- as.numeric(dir.new)
          lag <- as.numeric(lag.new)
        } else {
          if (ks.new == ks && value1 > value2) {
            value2 <- value1
            dir <- as.numeric(dir.new)
            lag <- as.numeric(lag.new)
          } #if
        } #else
      } #if
    } #for
  } else { #Value already achieved with k*s lags in the beginning
    value2 <- value
  }
  res <- list(lag_rect = lag, dir_rect = dir,
              lag_alllag = nrow(L), dir_alllag = alllag,
              lag_alldir = alldir, dir_alldir = ncol(L))
  res
}

# Lag-direction combinations with biggest values
lagdir_big <- function(L, thres = 0.8) {
  g <- 1
  value <- 0
  sortDT <- sort(L, decreasing = T)
  indx <- numeric(0)
  while (value < thres) {
    value <- sum(sortDT[1:g])
    indx <- cbind(indx, mtrindx(L, sortDT[g]))
    g <- g + 1    
  }
  index <- matrix(indx, ncol = 2, byrow = T)
  colnames(index) <- c("Lag", "Dir.")
  index
}

# Finding the matrix indices for lagdir_big
mtrindx <- function(i, val) {
  ind <- which(data.frame(i) == val)
  c.num <- ceiling(ind/nrow(i))
  r.num <- ind - nrow(i) * (c.num - 1)
  c(r.num, c.num)
}

