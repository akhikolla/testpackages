# Method AMUSEladle
AMUSEladle <- function(X, ...) UseMethod("AMUSEladle")

AMUSEladle.default <- function(X, tau = 1, l = 20, sim = c("geom", "fixed"), n.boot = 200,
                       ncomp = ifelse(ncol(X) > 10, floor(ncol(X)/log(ncol(X))), ncol(X) - 1), ...) {
  if (!is.numeric(X)) stop("non-numeric data")
  if (any(is.na(X) | is.infinite(X))) stop("missing/infinite values are not allowed")
  sim <- match.arg(sim)
  data.name <- deparse(substitute(X))
  method <- "AMUSE"

  n <- nrow(X)
  p <- ncol(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y 
  
  Mdata <- crossprod(Y[1:(n - tau), ], Y[(tau + 1):n, ])/(n - tau)
  Mdata.sym <- (Mdata + t(Mdata))/2
  Mdata.sym <- crossprod(Mdata.sym)

  EV.Mdata <- .Call("EIGEN", Mdata.sym, PACKAGE = "tsBSS")
  EVdata <- EV.Mdata$vectors
  EV.Mdatavalues <- as.vector(EV.Mdata$values)
  
  RES <- boot::tsboot(X, AMUSEbootLADLE, R = n.boot, sim = sim, l = l, EVdata = EVdata,
                tau = tau, rank = ncomp, ...)
  fis <- RES$t
  fn0 <- c(0, colMeans(fis))
  fn <- fn0/(1 + sum(fn0))
  phin <- EV.Mdatavalues[1:(ncomp + 1)]/(1 + sum(EV.Mdatavalues[1:(ncomp + 1)]))
  gn <- fn + phin
  est.k <- which.min(gn) - 1

  W <- crossprod(EVdata, prep$COV.sqrt.i)
  S <- ts(tcrossprod(prep$X.C, W))
  colnames(S) <- paste0("Series", 1:p)
    
  RES <- list(method = method, k = est.k, fn = fn, phin = phin, data.name = data.name,
              gn = gn, lambda = sort(EV.Mdata$values, decreasing = TRUE)[1:(ncomp + 1)],
              W = W, S = S, MU = prep$MEAN, sim = sim, lag = tau)
    
  class(RES) <- "ladle"
  RES
}


AMUSEladle.ts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- AMUSEladle.default(x, ...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}

AMUSEladle.xts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- AMUSEladle.default(x, ...)
  S <- xts::as.xts(RES$S)
  attr(S, "index") <- attr(X, "index")
  xts::xtsAttributes(S) <- xts::xtsAttributes(X) #attributes additional to zoo
  RES$S <- S
  RES
}

AMUSEladle.zoo <- function(X, ...) {
  x <- as.matrix(X)
  RES <- AMUSEladle.default(x, ...)
  S <- zoo::as.zoo(RES$S)
  attr(S, "index") <- attr(X, "index")
  RES$S <- S
  RES
}
