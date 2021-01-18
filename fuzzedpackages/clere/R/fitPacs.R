#' fitPacs function
#' 
#' This function implements the PACS (Pairwise Absolute Clustering and
#' Sparsity) methodology of Sharma DB et al. (2013). This methodology proposes
#' to estimate the regression coefficients by solving a penalized least squares
#' problem.  It imposes a constraint on Beta (the vector of regression
#' coefficients) that is a weighted combination of the L1 norm and the pairwise
#' L-infinity norm.  Upper-bounding the pairwise L-infinity norm enforces the
#' covariates to have close coefficients.  When the constraint is strong
#' enough, closeness translates into equality achieving thus a grouping
#' property.  For PACS, no software was available.  Only an R script was
#' released on Bondell's webpage
#' (http://www4.stat.ncsu.edu/~bondell/Software/PACS/PACS.R.r).  Since this R
#' script was running very slowly, we decided to reimplement it in C++ and
#' interfaced it with the present R package clere. This corresponds to the
#' option \code{type=1} in Bondell's script.
#' 
#' 
#' @param Y [numeric]: The vector of observed responses - size \code{n}.
#' @param X [matrix]: The matrix of predictors - size \code{n} rows and
#' \code{p} columns.
#' @param lambda [numeric]: A non-negative penalty term that controls
#' simultaneouly clusetering and sparsity.
#' @param betaInput [numeric]: A vector of initial guess of the model
#' parameters. The authors suggest to use coefficients obtained after fitting a
#' ridge regression with the shrinkage parameter selected using AIC criterion.
#' @param epsPACS [numeric]: A tolerance threshold that control the convergence
#' of the algorithm. The default value fixed in Bondell's initial script is
#' 1e-5.
#' @param nItMax [numeric]: Maximum number of iterations in the algorithm.
#' 
#' @return Object of class \code{\linkS4class{Pacs}} containing all the input
#' parameters plus parameter \code{a0} the intercept and parameter \code{K} the
#' dimensionality of the model.
#' 
#' @export
#' 
#' @seealso Overview : \code{\link{clere-package}} \cr 
#' Classes : \code{\linkS4class{Clere}}, \code{\linkS4class{Pacs}} \cr 
#' Methods : \code{\link{plot}}, \code{\link{clusters}}, \code{\link{predict}}, \code{\link{summary}} \cr 
#' Functions : \code{\link{fitClere}}, \code{\link{fitPacs}} 
#' Datasets : \code{\link{numExpRealData}}, \code{\link{numExpSimData}}, \code{\link{algoComp}}
#' 
#' @examples
#' 
#'    n     <- 100
#'    p     <-  20
#'    Beta  <- rep(c(0,2),10)
#'    eps   <- rnorm(n,sd=3)
#'    x     <- matrix(rnorm(n*p), nrow = n, ncol = p)
#'    y     <- as.numeric(10+x%*%Beta+eps)
#'    bInit <- lm(y~scale(x))$coefficients[-1]
#'    mod   <- fitPacs(Y=y,X=x,lambda=1.25,betaInput=bInit,epsPACS=1e-5,nItMax=1000) 
#' 
fitPacs <- function(
  Y,
  X,
  lambda = 0.5,
  betaInput,
  epsPACS = 1e-5,
  nItMax = 1000
) {
  if (lambda <= 0) {
    stop("[Clere:fitPacs] Non negative (or  = 0) values for tuning parameter lambda are not allowed!\n", call. = FALSE)
  }
  if (epsPACS <= 0) {
    stop("[Clere:fitPacs] Non negative (or  = 0) values for tolerance parameter epsPACS are not allowed!\n", call. = FALSE)
  }

  n <- nrow(X)
  p <- ncol(X)
  x <- scale(X)
  y <- Y - mean(Y)

  pacsObj <- methods::new("Pacs",
    y = y, x = x,
    n = as.integer(n), p = as.integer(p), nItMax = as.integer(nItMax),
    lambda = lambda, epsPACS = epsPACS, betaInput = betaInput
  )
  .Call("pacs", pacsObj, PACKAGE = "clere")
  littleeps <- 1e-7
  nround <- round(-log10(littleeps))
  rB <- round(pacsObj@betaOutput, nround)
  if (sum(rB == 0, na.rm = TRUE) > 0) {
    K <- length(unique(abs(rB[which(rB != 0)])))
  } else {
    K <- length(unique(abs(rB)))
  }
  pacsObj@a0 <- mean(Y - X %*% pacsObj@betaOutput)
  pacsObj@K <- K

  pacsObj
}
