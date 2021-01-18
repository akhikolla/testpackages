#---------------------------------------------
#' Generic S3 method vcov
#' @aliases vcov
#' @export
#' @param spbp a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the variance-covariance matrix associated the regression coefficients.
#'
vcov <- function(spbp, ...) UseMethod("vcov")

#---------------------------------------------
#' Covariance of the regression coefficients
#'
#' @aliases vcov.spbp
#' @rdname vcov-methods
#' @method vcov spbp
#' @export
#' @export vcov
#' @param spbp an object of the class spbp
#' @param ... further arguments passed to or from other methods.
#' @return  the variance-covariance matrix associated with the regression coefficients.
#'
#'
vcov.spbp <- function(spbp, ...){
  if(spbp$call$approach == "mle"){return(spbp$var)}
  else
    warning("not applicable, change approach to 'mle' to get covariance matrix")
}

#---------------------------------------------
#' Generic S3 method coef
#' @aliases coef
#' @export
#' @param spbp a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the estimated regression coefficients
#'
coef <- function(spbp, ...) UseMethod("coef")

#---------------------------------------------
#' Estimated regression coefficients
#'
#' @aliases coef.spbp
#' @rdname coef-methods
#' @method coef spbp
#' @export
#' @export coef
#' @param spbp an object of the class spbp
#' @param ... further arguments passed to or from other methods
#' @return  the estimated regression coefficients
#'
#'
coef.spbp <- function(spbp, ...){
  if(spbp$call$approach == "mle"){return(spbp$coefficients)}
  else
    warning("not applicable, change approach to 'mle' to get point estimates")
}


#---------------------------------------------
#' Generic S3 method confint
#' @aliases confint
#' @export
#' @param spbp a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the estimated regression coefficients
#'
confint <- function(spbp, ...) UseMethod("confint")

#---------------------------------------------
#' Confidence intervals for the regression coefficients
#'
#' @aliases confint.spbp
#' @rdname confint-methods
#' @method confint spbp
#' @export
#' @export confint
#' @param spbp an object of the class spbp
#' @param level the confidence level required
#' @param ... further arguments passed to or from other methods
#' @return  100(1-alpha) confidence intervals for the regression coefficients
#'
#'
confint.spbp <- function(spbp, level=0.95, ...){
  if(spbp$call$approach == "mle"){
    se <- sqrt(diag(vcov(spbp)))
    alpha <- 1-level
    CI <-   coef(spbp) +  se %o% c(-qnorm(1 - alpha/2),qnorm(1 - alpha/2))
    labels <- round(100*(c(alpha/2, 1-alpha/2)),1)
    colnames(CI) <- paste0(labels, "%")
    return(CI)
  }
  else
    warning("not applicable, change approach to 'mle' to get confidence intervals")
}


#' Bernstein basis polynomials calculations
#' @export
#' @param time a vector of times.
#' @param degree Bernstein polynomial degree
#' @param tau must be greater than times maximum value observed.
#' @return A list containing matrices b and B corresponding BP basis and corresponding tau value used to compute them.

bp.basis <- function(time, degree,  tau = max(time)){
  n <- length(time)

  ## error handling
  if(sum(time >= 0) != n)
    stop("time must be a positive vector.")

  if(degree < 0)
    stop("polynomial degree must be positive.")

  if(!(degree %% 1 == 0))
    stop("polynomial degree must be integer.")

  if(tau <  max(time))
    stop("tau must be greater than the last time.")

  k <- 1:degree
  b <- matrix(NA, n, degree )
  B <- matrix(NA, n, degree )
  y <- time/tau

  b <- sapply(k, function(k){dbeta(y, k, degree - k + 1) / tau})
  B <- sapply(k, function(k) pbeta(y, k, degree - k + 1) )

  # Equivalent to
  # for (i in 1:n){
  #   for(k in 1:degree){
  #     b[i,k] <- dbeta(y[i], k, degree - k + 1) / tau
  #     B[i,k] <- pbeta(y[i], k, degree - k + 1)
  #   }
  # }
  return(list(b = b, B = B, degree = degree, tau = tau))
}

#' Calculate the posterior mode
#' @export
#' @param ext rstan extracted sample.
#' @return A vector containing the posterior mode of each sample.

mode <- function(ext){
  f <- density(ext)
  pmode <- f$x[which.max(f$y)]
  return(pmode)
}

terms.inner <- function (x)
{
  if (class(x) == "formula")
    c(terms.inner(x[[2]]), terms.inner(x[[3]]))
  else if (class(x) == "call" && (x[[1]] != as.name("$") &&
                                  x[[1]] != as.name("["))) {
    if (x[[1]] == "+" || x[[1]] == "*" || x[[1]] == "-") {
      c(terms.inner(x[[2]]), terms.inner(x[[3]]))
    }
    else if (x[[1]] == as.name("Surv") || x[[1]] == as.name("rand"))
      unlist(lapply(x[-1], terms.inner))
    else terms.inner(x[[2]])
  }
  else (deparse(x))
}
####

blockSolve <- function(M, q){

  if(!is.matrix(M))
    stop("M is not a matrix")

  if(!is.numeric(q) || length(q) > 1 || q%%1 != 0)
    stop("q must be an integer")

  if(nrow(M) != ncol(M))
    stop("non square matrix")
  if(q == 0) return(M)

  n <- nrow(M)
  r <- (q+1)

  A <- M[1:q, 1:q];
  B <- M[1:q, r:n];
  C <- M[r:n,1:q];
  D <- M[r:n, r:n];

  S <- matrix(NA, nrow = nrow(M), ncol = ncol(M))

  ## aux op
  invA <- MASS::ginv(A)
  invPsis <-  MASS::ginv(D - C %*% invA %*% B)

  ##S11
  S[1:q, 1:q] <- invA + (invA %*% B %*% invPsis %*% C %*% invA)
  ##S12
  S[1:q, r:n] <- -(invA %*% B %*% invPsis)
  ##S21
  S[r:n,1:q] <- -(invPsis %*% C %*% invA)
  ##S22
  S[r:n, r:n] <- invPsis

  return(S)
}

solveAny <- function(A){
  tol <- .Machine$double.eps ## solve default tolerance
  class(S) <- "try-error"
  while(class(S) == "try-error"){
    S <- try(qr.solve(A, tol = tol), silent = T)
    tol <- tol^(1.00001)
  }
  return(S)
}

read_prior <- function(prior){
  aux <- unlist(strsplit(prior, "\\("))
  dist <- aux[1]
  aux2 <- unlist(strsplit(aux[2], "\\)"))[1]
  val <- unlist(strsplit(aux2, "\\,"))
  return(c(dist, val))
}
