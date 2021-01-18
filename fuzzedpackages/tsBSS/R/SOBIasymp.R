# Method SOBIasymp
SOBIasymp <- function(X, ...) UseMethod("SOBIasymp")

SOBIasymp.default <- function (X, k, tau = 1:12, eps = 1e-06, maxiter = 200, ...) {
  if (!is.numeric(X)) stop("non-numeric data")
  if (any(is.na(X) | is.infinite(X))) stop("missing/infinite values are not allowed")
  DNAME <- deparse(substitute(X))
  if (length(tau) == 1) tau <- 1:tau
  ntau <- length(tau)

  n <- nrow(X)
  p <- ncol(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y 
  
  R <- array(0, dim = c(p, p, ntau))
  n <- nrow(X)
  for (i in 1:ntau) {
    Yt <- Y[1:(n - tau[i]), ]
    Yti <- Y[(1 + tau[i]):n, ]
    Ri <- crossprod(Yt, Yti)/nrow(Yt)
    R[, , i] <- (Ri + t(Ri))/2
  }
  
  JD <- JADE::frjd(R, eps = eps, maxiter = maxiter)
  W <- crossprod(JD$V, prep$COV.sqrt.i)
  sumEVs <- apply(JD$D^2,1:2,sum)
  
  ORDER <- order(diag(sumEVs), decreasing = TRUE)
  D <- diag(sumEVs)[ORDER]
  W <- W[ORDER, ]
  Z <- tcrossprod(prep$X.C, W)
  
  D2 <- sumEVs[ORDER, ORDER]
  Tk <- n * sum((D2[(k + 1):p, (k + 1):p]))
  
  Z <- ts(Z)
  colnames(Z) <- paste0("Series", 1:p)
  
  names(Tk) <- "T"
  
  PARAMETER <- ntau/2 * (p - k) * ((p - k) + 1)
  names(PARAMETER) <- c("df")
  PVAL <- 1 - pchisq(Tk, PARAMETER)
  METHOD <- c("SOBI test for white noise processes")
  ALTERNATIVE <- paste0("there are fewer than ", p - k, 
                        " white noise components")
  RES <- list(statistic = Tk, p.value = PVAL, parameter = PARAMETER, 
               method = METHOD, data.name = DNAME, alternative = ALTERNATIVE, 
               k = k, W = W, S = Z, D = D, MU = prep$MEAN)
  class(RES) <- c("ictest", "htest")
  RES
}

SOBIasymp.ts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SOBIasymp.default(x, ...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}

SOBIasymp.xts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SOBIasymp.default(x, ...)
  S <- xts::as.xts(RES$S)
  attr(S, "index") <- attr(X, "index")
  xts::xtsAttributes(S) <- xts::xtsAttributes(X) #attributes additional to zoo
  RES$S <- S
  RES
}

SOBIasymp.zoo <- function(X, ...) {
  x <- as.matrix(X)
  RES <- SOBIasymp.default(x, ...)
  S <- zoo::as.zoo(RES$S)
  attr(S, "index") <- attr(X, "index")
  RES$S <- S
  RES
}

