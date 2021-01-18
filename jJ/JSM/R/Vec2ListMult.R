
Vec2ListMult <- function (para, ncz, ncb) {
  
  gamma <- para[1 : ncb]
  phi <- if (ncz > 0) para[(ncb + 1) : (ncb + ncz)] else numeric(0) # PZH change (ncw>0) to (ncz>0)
  alpha <- para[ncb + ncz + 1]
  Ysigma <- para[ncz + ncb + 2]
  Bsigma <- para[ncz + ncb + 3]

  result <- list(gamma = gamma, phi = phi, alpha = alpha, Ysigma = Ysigma, Bsigma = Bsigma)
  return(result)
}

