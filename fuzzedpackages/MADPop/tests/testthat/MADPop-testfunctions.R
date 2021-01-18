# sample from dirichlet distribution
rdirichlet <- function(n, alpha) {
  K <- length(alpha)
  X <- matrix(rgamma(n*K, shape = alpha), n, K, byrow = TRUE)
  sweep(X, 1, rowSums(X), "/")
}

# dirichlet-multinomial distribution
# x and alpha are vectors of the same length.
# just like dmultinom, this density is not vectorized (i.e., doesn't accept matrix inputs)
ddirmulti <- function(x, alpha, log = FALSE) {
  sx <- sum(x)
  sa <- sum(alpha)
  ans <- lgamma(sx+1) - sum(lgamma(x+1))
  ans <- ans + lgamma(sa) - lgamma(sx+sa)
  ans <- ans + sum(lgamma(x+alpha)) - sum(lgamma(alpha))
  if(!log) ans <- exp(ans)
  ans
}

# dirichlet-multinomial logliklihood
# vectorized to accept a matrix alpha and drops constants in x
# output length is nrow(alpha)
ldirmulti <- function(alpha, x) {
  sx <- sum(x)
  sa <- rowSums(alpha)
  ans <- lgamma(sa) - lgamma(sx+sa)
  ans + rowSums(lgamma(sweep(alpha, 2, x, "+")) - lgamma(alpha))
}

extract.iter <- function(stanfit, ii) {
  lapply(extract(stanfit), function(arr) {
    lst <- apply(arr, 1, list)
    lst[[ii]][[1]]
  })
}
