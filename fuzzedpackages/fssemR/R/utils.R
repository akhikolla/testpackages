## Objective value of multiple regression
##' @title obj.multiRegression
##' @param Xs eQTL matrices
##' @param Ys gene expression matrices
##' @param fit regression fit result object
##' @param trans  if rows for sample, trans = TRUE, otherwise, trans = FALSE. Default FALSE
##' @return error squared norm of ||(I-B)Y - FX||_2^2
obj.multiRegression = function(Xs, Ys, fit, trans = F) {
  m = length(Ys)
  if (is.matrix(Xs)) {
    Xs = lapply(1:m, function(i) {
      Xs
    })
  }
  if (!trans) {
    Xs = lapply(1:m, function(i) {
      t(Xs[[i]])
    })
    Ys = lapply(1:m, function(i) {
      t(Ys[[i]])
    })
  }
  error = .Call("ObjMultiReg", Xs, Ys, fit, PACKAGE = "fssemR")
  error
}


## center Xs and Ys
##' @title proc.centerFSSEM
##' @param Xs eQTL matrices
##' @param Ys list of gene expression matrices
##' @return  centered Xs and Ys and mean vectors
proc.centerFSSEM = function(Xs, Ys) {
  mX = rowMeans(Xs)
  mY = lapply(Ys, rowMeans)
  center = function(m) {
    apply(m, 1, function(x){ x - mean(x) })
  }
  Xs = center(Xs)
  Ys = lapply(Ys, center)
  list(Xs = t(Xs),
       Ys = lapply(Ys, t),
       mX = mX, mY = mY)
}

## logliklihood of FSSEM
##' @title logLikFSSEM
##' @param Bs  Network matrices
##' @param Wl  Weights for lasso term
##' @param Wf  Weights for fused term
##' @param lambda Hyperparameter of lasso term
##' @param rho Hyperparameter of fused lasso term
##' @param sigma2 noise variance
##' @param Dets determinants of I-B matrices
##' @param n number of observations
##' @param p number of genes
##' @return objective value of FSSEM with specified hyper-paramters
logLikFSSEM = function(Bs, Wl, Wf, lambda, rho, sigma2, Dets, n, p) {
  K = length(Bs)
  loglik = 0
  l1 = 0
  lf = 0
  rho = min(rho, 1e12)
  lf = rho * Wf * abs(Bs[[2]] - Bs[[1]])
  for (k in 1:K) {
    l1 = l1 + lambda * Wl[[k]] * abs(Bs[[k]])
    loglik = loglik - n[k] / 2 * log(Dets[k]**2)
  }
  diag(l1) = 0
  diag(lf) = 0
  loglik + p * sum(n) / 2 * log(sigma2) + sum(l1) + sum(lf)
}


##' function generator function
##' @title cwiseGradient4FSSEM
##' @param n number of observations
##' @param c cofactor vector
##' @param Y Matrix of gene expression
##' @param R Residual matrix
##' @param Y2norm Column of YtY
##' @param sigma2 noise variance
##' @return function whose argument is column vector bi
cwiseGradient4FSSEM = function(n, c, Y, R, Y2norm, sigma2) {
  function(x) {
    n * t(c) + (Y2norm * x - tcrossprod(R, Y)) / sigma2
  }
}

## lipschitz computation
# cwiseLipschitzFSSEMv0 = function(n, o, g, s, c2, Deti, Y2norm, sigma2, p) {
#  oxg = tcrossprod(o, g)
#  Imo = diag(p - 1) - g %*% oxg
#  sog = s %*% oxg
#  c = 1 - s %*% tcrossprod(o, s)
#  zta = 1e-12
#  x = -1 * tcrossprod(chol2inv(Imo + diag(p - 1) * zta), sog)
#  Li = n * c2 / (crossprod(x, Imo %*% x) + 2 * sog %*% x + c) / (Deti + 1e-16) + Y2norm / sigma2
#  while(Li < 0) {
#    zta = zta * 10
#    x = -1 * tcrossprod(chol2inv(Imo + diag(p - 1) * zta), sog)
#    Li = n * c2 / (crossprod(x, Imo %*% x) + 2 * sog %*% x + c) / (Deti + 1e-16) + Y2norm / sigma2
#  }
#  Li
#}

## speed-up lipschitz computation
cwiseLipschitz4FSSEM = function(n, z, c2, b2, Y2norm, sigma2) {
  n * c2 / (abs(z) - sqrt(c2 * b2))**2 + Y2norm / sigma2
}

## speed-up lipschitz computation
## abs(z - c %*% b) >= abs(z) - abs(c %*% b) >= abs(z) - sqrt(c2 * b2)
cwiseLipschitz4FSSEM2B = function(n, z, c2, b2, Y2norm, sigma2, ImB, i) {
  if (abs(z) - sqrt(c2 * b2) > 0) {
    n * c2 / (abs(z) - sqrt(c2 * b2))**2 + Y2norm / sigma2
  } else {
    cwiseLipschitzFSSEMv0(n, c2, Y2norm, sigma2, ImB, i)[1]
    ## Inf
  }
}

cwiseLipschitzFSSEMv0 = function(n, c2, Y2norm, sigma2, ImB, i) {
  ## x^tWx + 2R^tx + C
  p = nrow(ImB)
  I = chol2inv(chol(crossprod(ImB[, -i])))
  D = det(crossprod(ImB[, -i]))
  W = diag(p - 1) - ImB[-i, -i] %*% tcrossprod(I, ImB[-i, -i])
  R = ImB[-i, -i] %*% tcrossprod(I, ImB[i, -i, drop = FALSE])
  C = 1 - ImB[i, -i, drop = FALSE] %*% tcrossprod(I, ImB[i, -i, drop = FALSE])
  e = 1e-6
  v = solve(W + diag(e, p - 1), -1 * R)
  L = n * c2 / (crossprod(v, I %*% v) + 2 * crossprod(R, v) + C) / (D + e) + Y2norm / sigma2
  if (L < 0) {
    L = Inf
  }
  L
}


## Proximal operator for FLSA
proxFLSA = function(u, w, r, lambda, rho, c) {
  D   = u[[1]] - u[[2]]
  eQ  = (abs(D) <= 2 * rho * r / c)
  Df  = 1 - eQ
  rho = min(rho, 1e16)
  x = list(
    (u[[1]] + u[[2]]) / 2 * eQ + (u[[1]] - sign(D) * rho * r / c) * Df,
    (u[[1]] + u[[2]]) / 2 * eQ + (u[[2]] + sign(D) * rho * r / c) * Df
  )
  lapply(1:2, function(k) {
    xe = pmax(x[[k]] - lambda * (w[[1]] + w[[2]]) / 2 / c, 0) + pmin(x[[k]] + lambda * (w[[1]] + w[[2]]) / 2 / c, 0)
    xd = pmax(x[[k]] - lambda * w[[k]] / c, 0) + pmin(x[[k]] + lambda * w[[k]] / c, 0)
    xe * eQ + xd * Df
  })
}

sigma2FSSEM = function(Xs, Ys, Bs, Fs, n, p, m) {
  error = 0
  for(k in 1:m) {
    error = error + (norm(Ys[[k]] - Bs[[k]] %*% Ys[[k]] - Fs[[k]] %*% Xs, "f"))**2
  }
  error / (sum(n) * p)
}

inert_opt = function(opts = c("continuous", "linear"), init = 0) {
  switch(
    opts,
    "continuous" = function(k) {
      init
    },
    "linear"  = function(k) {
      (k - 1) / (k + 2)
    }
  )
}

##' @title inverseB
##' @description inverse matrices of B network for adaptive FSSEM
##' @param Bs list of network matrices
##' @return list of inversed B matrices
##' @export
inverseB = function(Bs) {
  lapply(Bs, function(B) { 1 / abs(B) })
}

##' @title invoneB
##' @description if you do not want to get inversed B matrces, invoneB gives you a matrix with constant 1 instead in FSSEM
##' @param Bs list of network matrices
##' @return list of invone B matrices
##' @export
invoneB = function(Bs) {
  lapply(Bs, function(B) { 1 / abs(sign(B)) })
}

##' @title flinvB
##' @description inversed difference of two B matrices. For adaptive fused lasso penalty
##' @param Bs list of network matrices
##' @return inversed difference matrices
##' @export
flinvB = function(Bs) {
  1 / abs(Bs[[2]] - Bs[[1]])
}

##' @title floneB
##' @description if you do not want adaptive fused lasso penalty, floneB replace flinvB
##' @param Bs list of network matrices
##' @return matrix whose entries are all 1
##' @export
floneB = function(Bs) {
  1 / abs(sign(Bs[[2]] - Bs[[1]]))
}

logLiklihood = function(Xs, Ys, Bs, Fs, mu, Dets, sigma2, p) {
  err  = 0
  logl = 0
  m = length(Ys)
  n = sapply(Ys, ncol)
  df = 0
  for(i in 1:m) {
    err  = err + norm(Ys[[i]] - Bs[[i]] %*% Ys[[i]] - Fs[[i]] %*% Xs - tcrossprod(mu[[i]], rep(1, n[i])), "f")**2
    logl = logl - n[i] / 2 * log(Dets[i]**2)
    df = df + sum(Bs[[2]] - Bs[[1]] != 0 & Bs[[i]] != 0) + p
  }
  df = df + sum(Bs[[2]] - Bs[[1]] == 0 & Bs[[1]] != 0) + 1
  sigma2 = err / (sum(n) * p)
  2 * (logl + sum(n) * p / 2 * log(2 * pi * sigma2)) + df * (log(n[1]) + log(n[2]))
}

bayesianInfocriterion = function(Xs, Ys, Bs, Fs, mu, Dets, sigma2, p) {
  logl = 0
  m = length(Ys)
  n = sapply(Ys, ncol)
  df = 0
  for(i in 1:m) {
    err  = norm(Ys[[i]] - Bs[[i]] %*% Ys[[i]] - Fs[[i]] %*% Xs - tcrossprod(mu[[i]], rep(1, n[i])), "f")**2
    logl = logl - n[i] / 2 * log(Dets[i]**2) + n[i] * p / 2 * log(sigma2) + err / (2 * sigma2)
    df = df + sum(Bs[[2]] - Bs[[1]] != 0 & Bs[[i]] != 0) + p
  }
  df = df + sum(Bs[[2]] - Bs[[1]] == 0 & Bs[[1]] != 0) + 1
  2 * logl + df * (log(n[1]) + log(n[2]))
}

L2lamax = function(Xs, Ys, Sk, n, p, k) {
  m = length(Ys)
  Xs = lapply(1:m, function(i) {
    t(Xs[[i]])
  })
  Ys = lapply(1:m, function(i) {
    t(Ys[[i]])
  })
  lambda = .Call("L2lamax", Xs, Ys, Sk, n, p, k, PACKAGE = "fssemR")
  lambda
}

##' @title TPR
##' @description Power of detection for network prediction
##' @param X list of predicted network matrices
##' @param B list of true network matrices
##' @param PREC precision threshold for FDR test. Default 0.
##' @export
TPR = function(X, B, PREC = 0) {
  X = as.matrix(X)
  B = as.matrix(B)
  sum(abs(X[B != 0]) > PREC) / sum(B != 0)
}

##' @title FDR
##' @description  False discovery rate for network prediction
##' @param X list of predicted network matrices
##' @param B list of true network matrices
##' @param PREC precision threshold for FDR test. Default 0.
##' @export
FDR = function(X, B, PREC = 0) {
  X = as.matrix(X)
  B = as.matrix(B)
  if(sum(X != 0) != 0) {
    sum(abs(X[B == 0]) > PREC) / sum(X != 0)
  } else {
    0
  }
}


## center Xs and Ys
##' @title proc.centerFSSEM2
##' @param Xs list of eQTL matrices
##' @param Ys list of gene expression matrices
##' @return  centered Xs and Ys and mean vectors
proc.centerFSSEM2 = function(Xs, Ys) {
  mX = lapply(Xs, rowMeans)
  mY = lapply(Ys, rowMeans)
  center = function(m) {
    apply(m, 1, function(x){ x - mean(x) })
  }
  Xs = lapply(Xs, center)
  Ys = lapply(Ys, center)
  list(Xs = lapply(Xs, t),
       Ys = lapply(Ys, t),
       mX = mX, mY = mY)
}


sigma2FSSEM2 = function(Xs, Ys, Bs, Fs, n, p, m) {
  error = 0
  for(k in 1:m) {
    error = error + (norm(Ys[[k]] - Bs[[k]] %*% Ys[[k]] - Fs[[k]] %*% Xs[[k]], "f"))**2
  }
  error / (sum(n) * p)
}

logLiklihood2 = function(Xs, Ys, Bs, Fs, mu, Dets, sigma2, p) {
  err  = 0
  logl = 0
  m = length(Ys)
  n = sapply(Ys, ncol)
  df = 0
  for(i in 1:m) {
    err  = err + norm(Ys[[i]] - Bs[[i]] %*% Ys[[i]] - Fs[[i]] %*% Xs[[i]] - tcrossprod(mu[[i]], rep(1, n[i])), "f")**2
    logl = logl - n[i] / 2 * log(Dets[i]**2)
    df = df + sum(Bs[[2]] - Bs[[1]] != 0 & Bs[[i]] != 0) + p
  }
  df = df + sum(Bs[[2]] - Bs[[1]] == 0 & Bs[[1]] != 0) + 1
  sigma2 = err / (sum(n) * p)
  2 * (logl + sum(n) * p / 2 * log(2 * pi * sigma2)) + df * (log(n[1]) + log(n[2]))
}

bayesianInfocriterion2 = function(Xs, Ys, Bs, Fs, mu, Dets, sigma2, p) {
  logl = 0
  m = length(Ys)
  n = sapply(Ys, ncol)
  df = 0
  for(i in 1:m) {
    err  = norm(Ys[[i]] - Bs[[i]] %*% Ys[[i]] - Fs[[i]] %*% Xs[[i]] - tcrossprod(mu[[i]], rep(1, n[i])), "f")**2
    logl = logl - n[i] / 2 * log(Dets[i]**2) + n[i] * p / 2 * log(sigma2) + err / (2 * sigma2)
    df = df + sum(Bs[[2]] - Bs[[1]] != 0 & Bs[[i]] != 0) + p
  }
  df = df + sum(Bs[[2]] - Bs[[1]] == 0 & Bs[[1]] != 0) + 1
  2 * logl + df * (log(n[1]) + log(n[2]))
}




