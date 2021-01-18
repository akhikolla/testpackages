#' @export
ridgeReg = function(X, Y, sk, gamma) {
  res = .Call("RidgeReg", t(X), t(Y), sk, gamma, PACKAGE = "SparseSEM")
  res
}


#' @export
multiReg = function(X, Y, sk, gamma, n, p, k) {
  m = length(Ys)
  Xs = lapply(1:m, function(i) {
    t(X)
  })
  Ys = lapply(1:m, function(i) {
    t(Ys[[i]])
  })
  res = .Call("MultiReg", Xs, Ys, sk, gamma, n, p, k, PACKAGE = "SparseSEM")
  res
}



#' @export
objvalReg = function(X, Y, fit) {
  m = length(Ys)
  Xs = lapply(1:m, function(i) {
    t(X)
  })
  Ys = lapply(1:m, function(i) {
    t(Ys[[i]])
  })
  error = .Call("ObjMultiReg", Xs, Ys, fit, PACKAGE = "SparseSEM")
  error
}

#' @export
objvalReg = function(X, Y, fit) {
  m = length(Ys)
  Xs = lapply(1:m, function(i) {
    t(Xs[[i]])
  })
  Ys = lapply(1:m, function(i) {
    t(Ys[[i]])
  })
  error = .Call("ObjMultiReg", Xs, Ys, fit, PACKAGE = "SparseSEM")
  error
}

#' @title SMLReg
#' @export
#' @examples
#' gamma = cv.multiRegression(data$obs$X, data$obs[c("Y1", "Y2")], data$obs$sk, ngamma = 50, nfold = 5, N, Ng, Nk)
#' fit = multiRegression(data$obs$X, data$obs[c("Y1", "Y2")], data$obs$sk, gamma, N, Ng, Nk, trans = F)
#' Xs = data$obs$X
#' Ys = data$obs[c("Y1", "Y2")]
#' Sk = data$obs$sk
#' fitc = SMLReg(Xs, Ys, Sk, fit$Bs, fit$Fs, fit$sigma2, lambda = 200, rho = 100, n = c(N, N), p = Ng, k = Nk)
SMLReg = function(X, Ys, Sk, Bs, Fs, sigma2, lambda, rho, n, p, k) {
  m = length(Ys)
  Xs = lapply(1:m, function(i) {
    t(X)
  })
  Ys = lapply(1:m, function(i) {
    t(Ys[[i]])
  })
  options = list(maxit = 100, threshold = 1e-4)
  fit = .Call("MultiFSSEMiPALM", Xs, Ys, Bs, Fs, Sk, sigma2, lambda, rho, options, n, p, k, m, PACKAGE = "SparseSEM")
  fit
}


#' @export
Lamax = function(X, Ys, Sk, n, p, k) {
  m = length(Ys)
  Xs = lapply(1:m, function(i) {
    t(X)
  })
  Ys = lapply(1:m, function(i) {
    t(Ys[[i]])
  })
  lambda = .Call("L2lamax", Xs, Ys, Sk, n, p, k, PACKAGE = "SparseSEM")
  lambda
}
