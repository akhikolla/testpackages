# The functions that are not used directly in package, but rather for testing C code.


ChangePointsId_R <- function(x) {
  len <- length(x)
  id <- seq_along(x)
  s <- sign(diff(x))
  id <- id[which(s != 0)]
  s <- diff(s[id])
  c(1, id[which(s != 0) + 1], len)
}


MovingSum <- function(x, n = 2) filter(x, rep(1, n), sides = 1)[-seq_len(n - 1)]


CheckSmallIntervalsB <- function(x, p, dn) {
  
  KeepB <- rep(TRUE, length(x))
  oKeepB <- NULL
  d <- 3
  while (!identical(KeepB, oKeepB) & d <= dn & sum(KeepB) > 3) {
    oKeepB <- KeepB
    d <- 3
    DropInt <- which(MovingSum((abs(diff(x[KeepB], lag = 1)))^p, d) < (abs(diff(x[KeepB], lag = d)))^p)
    if (length(DropInt) > 0) {
      Drop <- rep(DropInt, each = d - 1) + 1:(d - 1)
      KeepB[KeepB][Drop] <- FALSE
    }
    while (identical(KeepB, oKeepB) & d < dn & sum(KeepB) > d + 2) {
      d <- d + 2
      DropInt <- which(MovingSum((abs(diff(x[KeepB], lag = 1)))^p, d) < (abs(diff(x[KeepB], lag = d)))^p)
      if (length(DropInt) > 0) {
        Drop <- rep(DropInt, each = d - 1) + 1:(d - 1)
        KeepB[KeepB][Drop] <- FALSE
      }
    }
  }
  return(KeepB)
}

CheckSmallIntervals_R <- function(x, p, d) {
  id <- ChangePointsId_R(x)
  xx <- x[id]
  id[CheckSmallIntervalsB(xx, p, d)]
} 
