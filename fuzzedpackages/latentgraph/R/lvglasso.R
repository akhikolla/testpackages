lvglasso <- function(data, n, p, lambda1, lambda2, rule = "AND"){
  if (class(data)[1] == "matrix") {
    dat <- data
  }

  if (class(data)[1] == "data.frame") {
    dat <- as.matrix(data)
  }

  if (class(data)[1] == "list") {
    dat <- data[[1]]
    for (i in 2:n) {
      dat <- rbind(dat, data[[i]])
    }
  }

  Sigmahat <- var(dat)
  omega <- lvglasso_calc(Sigmahat,lambda1,lambda2,convergence=1e-10,maxiter=1000,rho=2.5)$S
  res <- omega
  # AND rule
  if (rule == "AND") {
    res <- res!=0
    res <- res*1
    res <- res + t(res)
    res[which(res==1)] <- 0
    res <- res/2
  }
  # OR rule
  if (rule == "OR") {
    res <- (res+t(res))/2
  }
  res <- res!=0
  res <- res*1

  return(list(omega=omega, theta = res, penalties = c(lambda1, lambda2)))
}
