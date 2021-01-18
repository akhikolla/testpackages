library(glmnet)
## yi[1 x n] = B_{i,-i}][1 x p-1]Yi[p-1 x n] + F[1 x sk]Xi[sk x n]
multiRegression2 = function(Xs, Ys, Sk, gamma, n, p, k, trans = FALSE, backend = "glmnet") {
  data = fssemR:::proc.centerFSSEM2(Xs, Ys)
  Xs = data$Xs
  Ys = data$Ys
  mX = data$mX
  mY = data$mY
  m = length(Ys)
  B = list(matrix(0, nrow = p, ncol = p), matrix(0, nrow = p, ncol = p))
  F = list(matrix(0, nrow = p, ncol = k), matrix(0, nrow = p, ncol = k))
  f = vector("list", 2)
  mu = vector("list", 2)
  err = 0
  for (j in 1:m) {
    f[[j]] = list()
    for (i in 1:p) {
      Xi = Xs[[j]][Sk[[i]], , drop = FALSE]      ## sk x n
      P  = diag(n) - crossprod(Xi, solve(tcrossprod(Xi))) %*% Xi  ## n x n
      yi = Ys[[j]][i, ,drop = FALSE]
      Yi = Ys[[j]][-i,]
      bi = glmnet(t(Yi %*% P), yi %*% P, alpha = 0, lambda = gamma / n, intercept = FALSE, standardize = FALSE)[["beta"]][,1]
      B[[j]][i, -i] = bi
      f[[j]][[i]] = (yi - bi %*% Yi) %*% crossprod(Xi, solve(tcrossprod(Xi)))
      F[[j]][i, Sk[[i]]] = f[[j]][[i]]
      err = err + tcrossprod(yi - bi %*% Yi - f[[j]][[i]] %*% Xi)
      f[[j]][[i]] = t(f[[j]][[i]])
    }
    mu[[j]] = (diag(p) - B[[j]]) %*% mY[[j]] - F[[j]] %*% mX[[j]]
    F[[j]] = t(F[[j]])
  }
  sigma2 = err / 2 / (n * p - 1)
  list(
    Bs = B,
    fs = f,
    Fs = F,
    mu = mu,
    sigma2 = sigma2
  )
}
