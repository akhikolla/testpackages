# Method gSOBI
gSOBI <- function(X, ...) UseMethod("gSOBI")

# gSOBI function for PVC (combination of SOBI and vSOBI)
gSOBI.default <- function(X, k1 = 1:12, k2 = 1:3, b = 0.9, eps = 1e-06, maxiter = 1000, 
                   ordered = FALSE, acfk = NULL, original = TRUE, alpha = 0.05, ...) {
  if (!is.numeric(X)) stop("non-numeric data")
  if (any(is.na(X) | is.infinite(X))) stop("missing/infinite values are not allowed")
  n <- nrow(X)
  p <- ncol(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y 
  
  U <- diag(p) #Initial value for the orthogonal matrix U
  crit <- Inf
  iter <- 0
  K1 <- length(k1)
  K2 <- length(k2)
  Tk1 <- array(NA, dim = c(p, p, K1))
  Tk2 <- array(NA, dim = c(p, p, K2))
  while (crit > eps) {
    for (i in 1:K1) {
      Tk1[ , , i] <- .Call("TIK1", Y, U, k1[i], PACKAGE = "tsBSS")$Tik
    }
    for (i in 1:K2) {
      Tk2[ , , i] <- .Call( "TIK", Y, U, k2[i], method = 3, PACKAGE = "tsBSS")$Tik
    }
    TU <- b*apply(Tk1, c(1, 2), sum) + (1 - b)*apply(Tk2, c(1, 2), sum)
    EVDt <- .Call("EIGEN", tcrossprod(TU), PACKAGE = "tsBSS")
    COVt.sqrt.i <- EVDt$vectors %*% tcrossprod(diag(as.vector(EVDt$values)^(-0.5)), EVDt$vectors)
    U.new <- COVt.sqrt.i %*% TU #Updated U
    crit <- sqrt(sum((abs(U.new) - abs(U))^2)) #Comparing the current and the new matrix U
    iter <- iter + 1
    if (iter > maxiter) stop("maxiter reached without convergence")
    U <- U.new
  } #While the criterion value is below tolerance value
  W <- crossprod(U, prep$COV.sqrt.i) #Unmixing matrix
  S <- tcrossprod(prep$X.C, W)
  if (ordered == TRUE) { #Ordering by volatility
    if (is.null(acfk) == TRUE) { acfk <- 1:max(k1, k2) }
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
  RES <- list(W = W, k1 = k1, k2 = k2, S = S, MU = prep$MEAN)
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

gSOBI.ts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- gSOBI.default(x, ...)
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

gSOBI.xts <- function(X, ...) {
  x <- as.matrix(X)
  RES <- gSOBI.default(x, ...)
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

gSOBI.zoo <- function(X, ...) {
  x <- as.matrix(X)
  RES <- gSOBI.default(x, ...)
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

#Function to order the components (calculates the volatilities of the components or their residuals)
# Used by the functions gFOBI, gJADE, gSOBI, vSOBI, PVC and FixNA
ordf <- function(S, acfk, p, W, alpha, ...) {
  lblin1 <- lbtest(S, acfk, "linear") #Linear autocorrelations exist?
  lineff <- lblin1$p_val
  armaeff <- (lineff < alpha) #Logical vector: TRUE if series is replaced with residuals; FALSE if not

  S2 <- S
  fits <- vector("list", p)
  # Replaces series with residuals (if autocorrelation exists)
  fits[armaeff == 1] <- lapply(as.data.frame(S[, armaeff == 1]), forecast::auto.arima, stationary = TRUE, seasonal = FALSE, ...)
  S2[, armaeff == 1] <- sapply(fits[armaeff == 1], residuals)
  
  vol <- lbtest(S2, acfk, "squared")
  ord <- vol$TS
  ordered <- order(ord, decreasing = TRUE)
  
  list(S = S[, ordered],
       RS = S2[, ordered], #Residuals if ARMA effects; otherwise original independent component
       W = W[ordered, ],
       fits = fits[ordered],
       volTS = vol$TS[ordered], volP = vol$p_val[ordered],
       armaeff = armaeff[ordered], linTS = lblin1$TS[ordered], linP = lineff[ordered])
}

`print.bssvol` <- function(x, ...) {
  print.listof(x[(names(x) != "S") & (names(x) != "Sraw") & (names(x) != "MU") & (names(x) != "fits") & (names(x) != "armaeff") & (names(x) != "linTS")
                 & (names(x) != "linP") & (names(x) != "volTS") & (names(x) != "volP")], ...)
}

`plot.bssvol` <- function(x, ...) {
  S <- x$S
  if(ncol(S) <= 2) {
    plot(S, ...)
    } else {
      if (any(class(S) %in% c("mts", "xts", "zoo"))) {
        plot(S, ...)
        } else {
          pairs(S, ...)
        }
    }
}
