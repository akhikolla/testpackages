# Ljung-Box-type test
# Option "linear": Checking if linear autocorrelations are present
#  - whether ARMA-components are present or not
# Option "squared": Checking if second order autocorrelations are present
# - whether volatility is constant or not

lbtest <- function(X, k, type = c("squared", "linear")) {
  if (!is.numeric(X)) stop("non-numeric data")
  if (any(is.na(X) | is.infinite(X))) stop("missing/infinite values are not allowed")
  Xname <- deparse(substitute(X))
  type <- match.arg(type)
  if (is.null(ncol(X))) {
    X <- as.matrix(X) #if X is given as a vector
    colnames(X) <- "Series 1"
  }
  TS <- switch(type, 
                linear = as.vector(.Call( "lblinM", X, k, PACKAGE = "tsBSS")$RES),
                squared = as.vector(.Call("lbsqM", X, k, PACKAGE = "tsBSS")$RES))
  K <- length(k)
  p_val <- 1 - pchisq(TS, K)
  RES <- list(TS = TS, p_val = p_val, Xname = Xname,
              varnames = colnames(X), k = k, K = K, type = type)
  class(RES) <- "lbtest"
  RES
}

`print.lbtest` <- function(x, digits = 3, ...) {
  cat("\n")
  cat(strwrap(paste("Serial autocorrelation test for", x$Xname), prefix = "\t"), sep = "\n")
  cat("\n")
  cat(paste("Testing for", x$type), "autocorrelations based on lag(s) ")
  cat(x$k, "\n")
  cat(paste("Based on a chi squared test with", x$K), "degrees of freedom \n\n")
  cat("The test statistic and the corresponding p-value for each series: \n")
  # If no variable names available, names will be Series 1, Series 2, ...
  if (is.null(x$varnames)) x$varnames <- paste("Series", 1:length(x$TS))
  restab <- as.data.frame(cbind(x$varnames, round(x$TS, digits), format.pval(x$p_val, digits = digits)))
  colnames(restab) <- c("Series", "Statistic", "p-value")
  print(restab, row.names = FALSE, digits = 4)
  cat("\n")
  cat("For each series the alternative hypothesis is:", paste("serial", x$type), "correlation exists \n")
}



