# Method vSOBI
vSOBI <- function(X, ...) UseMethod("vSOBI")

# main function for vSOBI
vSOBI.default <- function(X, k = 1:12, eps = 1e-06, maxiter = 1000, G = c("pow", "lcosh"), 
                           ordered = FALSE, acfk = NULL, original = TRUE, alpha = 0.05, ...) {
  if (!is.numeric(X)) stop("non-numeric data")
  if (any(is.na(X) | is.infinite(X))) stop("missing/infinite values are not allowed")
  G <- match.arg(G)

  n <- nrow(X)
  p <- ncol(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y 
  
  U <- diag(p) #Initial value for the orthogonal matrix U
  crit <- Inf
  iter <- 0
  if (length(which(k < 1) != 0)) stop("only non-zero lags allowed")
  nk <- length(k)
  Tk <- array(NA, dim = c(p, p, nk))
  while (crit > eps) {
    if (G == "pow") {
      for (i in 1:nk) {
        Tk[ , , i] <- .Call( "TIK", Y, U, k[i], method = 3, PACKAGE = "tsBSS")$Tik
      }
    } else {
      for (i in 1:nk) {
        Tk[ , , i] <- .Call( "TIKlc", Y, U, k[i], method = 3, PACKAGE = "tsBSS")$Tik
      }
    }
    TU <- apply(Tk, c(1, 2), sum)
    EVDt <- .Call("EIGEN", tcrossprod(TU), PACKAGE = "tsBSS")
    COVt.sqrt.i <- EVDt$vectors %*% tcrossprod(diag(as.vector(EVDt$values)^(-0.5)), EVDt$vectors)
    U.new <- COVt.sqrt.i %*% TU #Updated U
    crit <- sqrt(sum((abs(U.new) - abs(U))^2)) #Comparing the current and the new matrix U
    iter <- iter + 1
    if (iter > maxiter) stop("maxiter reached without convergence")
    U <- U.new
  } #While the criterion value is below tolerance value.
  
  W <- crossprod(U, prep$COV.sqrt.i) #Unmixing matrix
  S <- tcrossprod(prep$X.C, W)
  if (ordered == TRUE) { #Ordering by volatility
    if (is.null(acfk) == TRUE) { acfk <- k }
    ord <- ordf(S, acfk, p, W, alpha, ...)
    W <- ord$W
    if (original == TRUE) {
      S <- ord$S # Original independent components
    } else {
      S <- ord$RS # Residuals based on ARMA fit, if applicable; otherwise original IC's
      Sraw <- ord$S
      Sraw <- ts(Sraw, names = paste("Series", 1:p))
    }
  }
  S <- ts(S, names = paste("Series", 1:p))
  RES <- list(W = W, k = k, S = S, MU = prep$MEAN)
  if (ordered == TRUE) {
    if (original == FALSE) {
      RES$Sraw <- Sraw
    }
    RES$fits <- ord$fits
    RES$armaeff <- ord$armaeff
    RES$linTS <- ord$linTS
    RES$linP <- ord$linP
    RES$volTS <- ord$volTS
    RES$volP <- ord$volP
  }
  class(RES) <- c("bssvol", "bss")
  RES
}


vSOBI.ts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- vSOBI.default(x, ...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  if (!is.null(RES$Sraw)) {
    Sraw <- RES$Sraw
    attr(Sraw, "tsp") <- attr(X, "tsp")
    RES$Sraw <- Sraw
  }
  RES
}

vSOBI.xts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- vSOBI.default(x, ...)
  S <- xts::as.xts(RES$S)
  attr(S, "index") <- attr(X, "index")
  xts::xtsAttributes(S) <- xts::xtsAttributes(X) #attributes additional to zoo
  RES$S <- S
  if (!is.null(RES$Sraw)) {
    Sraw <- xts::as.xts(RES$Sraw)
    attr(Sraw, "index") <- attr(X, "index")
    xts::xtsAttributes(Sraw) <- xts::xtsAttributes(X)
    RES$Sraw <- Sraw
  }
  RES
}

vSOBI.zoo <- function(X, ...) {
  x <- as.matrix(X)
  RES <- vSOBI.default(x, ...)
  S <- zoo::as.zoo(RES$S)
  attr(S, "index") <- attr(X, "index")
  RES$S <- S
  if (!is.null(RES$Sraw)) {
    Sraw <- zoo::as.zoo(RES$Sraw)
    attr(Sraw, "index") <- attr(X, "index")
    RES$Sraw <- Sraw
  }
  RES
}