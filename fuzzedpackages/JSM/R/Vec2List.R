
Vec2List <- function (para, ncx, ncz, ncw) {
  
  p <- ncz * (ncz + 1) / 2
  beta <- para[1 : ncx]
  phi <- if (ncw > 0) para[(ncx + 1) : (ncx + ncw)] else numeric(0)
  alpha <- para[ncx + ncw + 1]
  Ysigma <- para[ncx + ncw + 2]
  Bsigma <- para[(ncx + ncw + 3):(ncx + ncw + p + 2)]
  if (ncz > 1) {
    Bsigma.new <- matrix(0, ncz, ncz)
    Bsigma.new[lower.tri(Bsigma.new, diag = TRUE)] <- Bsigma
    Bsigma.new <- Bsigma.new + t(Bsigma.new) - diag(diag(Bsigma.new))
    Bsigma <- Bsigma.new
  }
  result <- list(beta = beta, phi = phi, alpha = alpha, Ysigma = Ysigma, Bsigma = Bsigma)
  return(result)
}

