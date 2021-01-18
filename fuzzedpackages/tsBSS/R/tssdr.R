#Time series supervised dimension reduction (tssdr)
# Method tssdr
tssdr <- function(y, X, ...) UseMethod("tssdr")


# main function for tssdr (time series supervised dimesion reduction)
tssdr.default <- function(y, X, algorithm = c("TSIR", "TSAVE", "TSSH"), k = 1:12, H = 10, weight = 0.5, method = c("frjd", "rjd"),
                   eps = 1e-06, maxiter = 1000, ...) {
  if (!is.numeric(X)) stop("non-numeric data")
  if (any(is.na(X) | is.infinite(X))) stop("missing/infinite values are not allowed")
  
  #weight only applies if algorithm = "TSSH"
  if (length(y) != nrow(X)) stop("y and X have different lengths!")
  algorithm <- match.arg(algorithm)
  method <- match.arg(method)
  if ((algorithm == "TSSH") & length(H) == 1) {
    H <- rep(H, 2)
    warning("H should be a 2-vector for TSSH. Using the given H for both parts.")
  }
  if ((algorithm != "TSSH") & length(H) == 2) {
    stop('H should be a scalar for TSIR and TSAVE!')
  }
  
  n <- nrow(X)
  p <- ncol(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y 
  
  nk <- length(k)
  if (length(H) == 1) { #For TSIR and TSAVE
    slices <- as.matrix(cut(y, breaks = c(quantile(y, probs = seq(0, 1, by = 1/H))),
                          include.lowest = TRUE, labels = FALSE))
  } else { #For TSSH
    slices <- as.matrix(cut(y, breaks = c(quantile(y, probs = seq(0, 1, by = 1/H[1]))),
                            include.lowest = TRUE, labels = FALSE))
    slices2 <- as.matrix(cut(y, breaks = c(quantile(y, probs = seq(0, 1, by = 1/H[2]))),
                            include.lowest = TRUE, labels = FALSE))
  }
  R <- array(0, dim = c(p, p, nk))
  switch(algorithm,
         TSIR = {
           for (i in 1:length(k)) {
             R[, , i] <- .Call("TSIR", Y, slices, k[i], H, PACKAGE = "tsBSS")$RES
           }
         },
         TSAVE = {
           for (i in 1:length(k)) {
             R[, , i] <- .Call("TSAVE", Y, slices, k[i], H, PACKAGE = "tsBSS")$RES
           }
         },
         TSSH = {
           for (i in 1:length(k)) {
             R[, , i] <- (1 - weight)*.Call("TSIR", Y, slices, k[i], H[1], PACKAGE = "tsBSS")$RES +
              weight*.Call("TSAVE", Y, slices2, k[i], H[2], PACKAGE = "tsBSS")$RES
           }
         }
  )
  # Joint diagonalization
  JD <- switch(method, frjd = {
    JADE::frjd(R, eps = eps, maxiter = maxiter, ...)
  }, rjd = {
    JADE::rjd(R, eps = eps, maxiter = maxiter, ...)
  })
  sumdg <- diag(apply(JD$D, 1:2, sum))
  ord <- order(sumdg, decreasing = TRUE)
  P <- diag(p)
  P <- P[ord, ]
  D <- JD$D
  DTable <- matrix(0, ncol = p, nrow = nk)
  for (j in 1:length(k)) {
    D[ , , j] <- P %*% tcrossprod(D[ , , j], P) #Diagonal elements are now in order
    DTable[j, ] <- diag(D[ , , j])
  }
  colnames(DTable) <- paste("Dir.", (1:p), sep = "")
  rownames(DTable) <- paste("Lag ", k, sep = "")
  V <- JD$V %*% t(P)
  W <- crossprod(V, prep$COV.sqrt.i) #p*p matrix
  S <- tcrossprod(prep$X.C, W) #Values for the all possible directions
  S <- ts(cbind(y, S)) # Response included
  colnames(S) <- c("y", paste("Series", 1:p))
  RES <- list(W = W, k = k, S = S, L = DTable/sum(DTable), H = H,
              yname = deparse(substitute(y)),
              Xname = deparse(substitute(X)),
              algorithm = algorithm)
  class(RES) <- "tssdr"
  RES
}

tssdr.ts <- function(y, X, ...) {
  yy <- as.vector(y)
  x <- as.matrix(X)
  RES <- tssdr.default(yy, x, ...)
  S <- RES$S
  attr(S, "tsp") <- attr(X, "tsp")
  RES$S <- S
  RES
}

tssdr.xts <- function(y, X, ...) {
  yy <- as.vector(y)
  x <- as.matrix(X)
  RES <- tssdr.default(yy, x, ...)
  S <- xts::as.xts(RES$S)
  attr(S, "index") <- attr(X, "index")
  xts::xtsAttributes(S) <- xts::xtsAttributes(X) #attributes additional to zoo
  RES$S <- S
  RES
}

tssdr.zoo <- function(y, X, ...) {
  yy <- as.vector(y)
  x <- as.matrix(X)
  RES <- tssdr.default(yy, x, ...)
  S <- zoo::as.zoo(RES$S)
  attr(S, "index") <- attr(X, "index")
  RES$S <- S
  RES
}

# Printing method for objects of class "tssdr"
`print.tssdr` <- function(x, digits = 3, ...) {
  print.listof(x[(names(x) == "W") | (names(x) == "k") | (names(x) == "L")], digits = digits, ...)
}

# Extracting the estimated directions:
# Components method for objects of class "tssdr"
`components.tssdr` <- function(x, ...) x$S[, -1]

# Plotting method for objects of class "tssdr" (R's basic time series plot)
`plot.tssdr` <- function(x, main = "The response and the directions", ...) {
  S <- x$S
  if(ncol(S) <= 2) {
    plot(S, main = main, ...)
  } else {
    if (any(class(S) %in% c("mts", "xts", "zoo"))) {
      plot(S, main = main, ...)
    } else {
      pairs(S, main = main, ...)
    }
  }
}