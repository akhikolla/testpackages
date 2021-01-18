
Indexing <- function (n) {
  
  p <- n * (n + 1) / 2
  mtx <- matrix(0, n, n)
  mtx[lower.tri(mtx, diag = TRUE)] <- 1:p
  mtx <- mtx + t(mtx) - diag(diag(mtx))
  return(c(mtx))
}
