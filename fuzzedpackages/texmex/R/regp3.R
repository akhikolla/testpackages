#' @rdname degp3
#' @export
regp3 <- function(n, kappa=1, sigma, xi, u=0){

  kappa <- rep(kappa, length.out=n)
  sigma <- rep(sigma, length.out=n)
  xi <- rep(xi, length.out=n)
  u <- rep(u, length.out=n)

  qegp3(-rexp(n), kappa, sigma, xi, u, log.p=TRUE)
}


