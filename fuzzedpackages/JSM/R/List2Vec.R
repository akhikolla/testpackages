
List2Vec <- function (theta) {

  Bsigma <- theta$Bsigma
  Bsigma <- if (is.matrix(Bsigma)) Bsigma[lower.tri(Bsigma, diag = TRUE)] else Bsigma
  
  para <- c(theta$beta, theta$phi, theta$alpha, theta$Ysigma, Bsigma)
  return(para)
}
