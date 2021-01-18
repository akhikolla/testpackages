# Method gJADE
gJADE <- function(X, ...) UseMethod("gJADE")

# main function for gJADE
gJADE.default <- function(X, k = 0:12, eps = 1e-06, maxiter = 100, method = c("frjd", "rjd"), 
                          na.action = na.fail, weight = NULL,
                          ordered = FALSE, acfk = NULL, original = TRUE, alpha = 0.05, ...) {
  if (!is.numeric(X)) stop("non-numeric data")
  if (any(is.na(X) | is.infinite(X))) stop("missing/infinite values are not allowed")
  nk <- length(k)
  method <- match.arg(method)

  n <- nrow(X)
  p <- ncol(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y 
  
  ccks <- .Call( "CCK", Y, k, PACKAGE = "tsBSS")$CCK
  U <- switch(method, frjd = {
    JADE::frjd(ccks, eps = eps, maxiter = maxiter, na.action = na.action, weight = weight)$V
  }, rjd = {
    JADE::rjd(ccks, eps = eps, maxiter = maxiter, na.action = na.action)$V
  })
  W <- crossprod(U, prep$COV.sqrt.i)
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

gJADE.ts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- gJADE.default(x, ...)
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

gJADE.xts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- gJADE.default(x, ...)
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

gJADE.zoo <- function(X, ...) {
  x <- as.matrix(X)
  RES <- gJADE.default(x, ...)
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