# Method SOBIladle
SOBIladle <- function(X, ...) UseMethod("SOBIladle")

SOBIladle.default <- function(X, tau = 1:12, l = 20, sim = c("geom", "fixed"), n.boot = 200,
                      ncomp = ifelse(ncol(X) > 10, floor(ncol(X)/log(ncol(X))), ncol(X) - 1),
                      maxiter = 1000, eps = 1e-06, ...) {
  if (!is.numeric(X)) stop("non-numeric data")
  if (any(is.na(X) | is.infinite(X))) stop("missing/infinite values are not allowed")
  sim <- match.arg(sim)
  data.name <- deparse(substitute(X))
  method <- "SOBI"
  if (length(tau) == 1) tau <- 1:tau
  
  n <- nrow(X)
  p <- ncol(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y
  
  M_array <- array(0, dim = c(p, p, length(tau)))
  for (t in 1:length(tau)) {
    M_array[, , t] <- crossprod(Y[1:(n - tau[t]), ], Y[(tau[t] + 1):n, ])/(n - tau[t])
    M_array[, , t] <- (M_array[, , t] + t(M_array[, , t]))/2
  }

  EV.Mdata <- MSobi(X, k_set = tau)
  frjddata <- JADE::frjd.int(EV.Mdata, eps = eps, maxiter = maxiter)
  Dfrjddata <- diag(apply(frjddata$D^2, 1:2, sum))
  EVdata <- frjddata$V[, order(Dfrjddata, decreasing = TRUE)]

  RES <- boot::tsboot(X, SOBIbootLADLE, R = n.boot, sim = sim, l = l, EVdata = EVdata, tau = tau,
                      rank = ncomp, maxiter = maxiter, eps = eps, ...)

  fis <- RES$t
  fn0 <- c(0, colMeans(fis))
  fn <- fn0/(1 + sum(fn0))
  phin <- sort(Dfrjddata, decreasing = TRUE)[1:(ncomp + 1)]/(1 + sum(sort(Dfrjddata, decreasing = TRUE)[1:(ncomp + 1)]))
  gn <- fn + phin
  est.k <- which.min(gn) - 1
  W <- crossprod(EVdata, prep$COV.sqrt.i)
  S <- tcrossprod(prep$X.C, W)
  S <- ts(tcrossprod(prep$X.C, W))
  colnames(S) <- paste0("Series", 1:p)

  RES <- list(method = method, k = est.k, fn = fn, phin = phin, data.name = data.name,
      gn = gn, lambda = sort(Dfrjddata, decreasing = TRUE)[1:(ncomp + 1)],
      W = W, S = S, MU = prep$MEAN, sim = sim, tau = tau)
  
  class(RES) <- "ladle"
  RES
  
}

SOBIladle.ts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SOBIladle.default(x, ...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}

SOBIladle.xts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SOBIladle.default(x, ...)
  S <- xts::as.xts(RES$S)
  attr(S, "index") <- attr(X, "index")
  xts::xtsAttributes(S) <- xts::xtsAttributes(X) #attributes additional to zoo
  RES$S <- S
  RES
}

SOBIladle.zoo <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SOBIladle.default(x, ...)
  S <- zoo::as.zoo(RES$S)
  attr(S, "index") <- attr(X, "index")
  RES$S <- S
  RES
}
