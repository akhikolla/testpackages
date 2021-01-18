# Method AMUSEasymp
AMUSEasymp <- function(X, ...) UseMethod("AMUSEasymp")

AMUSEasymp.default <- function (X, k, tau = 1, ...) {
  if (!is.numeric(X)) stop("non-numeric data")
  if (any(is.na(X) | is.infinite(X))) stop("missing/infinite values are not allowed")
  DNAME <- deparse(substitute(X))
  k <- abs(as.integer(k))
  tau <- abs(as.integer(tau))

  n <- nrow(X)
  p <- ncol(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y 
  
  Yt <- Y[1:(n - tau), ]
  Yti <- Y[(1 + tau):n, ]
  R <- crossprod(Yt, Yti)/nrow(Yt)
  R <- (R + t(R))/2
  
  EVD <- .Call("EIGEN", R, PACKAGE = "tsBSS")
  EVDvalues <- as.vector(EVD$values)
  W <- crossprod(EVD$vectors, prep$COV.sqrt.i)
  
  ORDER <- order(EVDvalues^2, decreasing = TRUE)
  D <- EVDvalues[ORDER]
  W <- W[ORDER, ]
  Z <- tcrossprod(prep$X.C, W)
  
  D2 <- D^2
  Tk <- n * sum((D2[(k + 1):p]))
  
  Z <- ts(Z)
  colnames(Z) <- paste0("Series", 1:p)
  
  names(Tk) <- "T"
  PARAMETER <- 0.5 * (p - k) * ((p - k) + 1)
  names(PARAMETER) <- c("df")
  PVAL <- 1 - pchisq(Tk, PARAMETER)
  METHOD <- c("SOBI test for white noise processes")
  ALTERNATIVE <- paste0("there are fewer than ", p - k, 
                        " white noise components")
  RES <- list(statistic = Tk, p.value = PVAL, parameter = PARAMETER, 
               method = METHOD, data.name = DNAME, alternative = ALTERNATIVE, 
               k = k, W = W, S = Z, D = D, MU = prep$MEAN, tau = tau)
  class(RES) <- c("ictest", "htest")
  RES
}

AMUSEasymp.ts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- AMUSEasymp.default(x, ...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}

AMUSEasymp.xts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- AMUSEasymp.default(x, ...)
  S <- xts::as.xts(RES$S)
  attr(S, "index") <- attr(X, "index")
  xts::xtsAttributes(S) <- xts::xtsAttributes(X) #attributes additional to zoo
  RES$S <- S
  RES
}

AMUSEasymp.zoo <- function(X, ...) {
  x <- as.matrix(X)
  RES <- AMUSEasymp.default(x, ...)
  S <- zoo::as.zoo(RES$S)
  attr(S, "index") <- attr(X, "index")
  RES$S <- S
  RES
}
