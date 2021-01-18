
List2VecMult <- function (theta) {
  para <- c(theta$gamma, theta$phi, theta$alpha, theta$Ysigma, theta$Bsigma)
  return(para)
}
