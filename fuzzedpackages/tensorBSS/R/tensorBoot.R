tensorBoot <- function(x, replace = TRUE) {
  dimx <- dim(x)
  n1 <- prod(dimx[-length(dimx)])
  n2 <- dimx[length(dimx)]
  dim(x) <- c(n1, n2)
  x <- x[, sample(1:n2, n2, replace = replace)]
  dim(x) <- dimx
  return(x)
}
