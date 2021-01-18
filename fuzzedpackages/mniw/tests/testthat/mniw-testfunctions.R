########################

# R functions for the mniw package

############################

# matrix of iid N(0,1)
rMnorm <- function(p, q) {
  if(missing(q)) q <- p
  matrix(rnorm(p*q), p, q)
}

# log-determinant
ldet <- function(X) {
  determinant(X, log = TRUE)$modulus[1]
}
ldetV <- function(V) {
  solveV(V, x = rep(1, nrow(V)), ldV = TRUE)$ldV
}

# log multi-gamma function
lmgamma <- function(x, p) p*(p-1)/4 * log(pi) + sum(lgamma(x + (1-1:p)/2))

# density of wishart and inverse wishart
dwishR <- function(X, Psi, nu, inverse = FALSE, log = FALSE) {
  if(length(X) == 1) {
    X <- as.matrix(X)
    Psi <- as.matrix(Psi)
  }
  d <- nrow(X)
  if(!inverse) {
    ans <- (nu-d-1) * ldetV(X)
    ans <- ans - nu * ldetV(Psi)
    ans <- ans - sum(diag(solveV(Psi, X)))
  } else {
    ans <- -(nu+d+1) * ldetV(X)
    ans <- ans + nu * ldetV(Psi)
    ans <- ans - sum(diag(solveV(X, Psi)))
  }
  ans <- ans - nu*d * log(2)
  ans <- 0.5 * ans - lmgamma(nu/2, d)
  if(!log) ans <- exp(ans)
  ans
}

# simulation of wishart/inverse wishart
rwishR <- function(Psi, nu, inverse = FALSE) {
  q <- nrow(Psi)
  XL <- matrix(0, q,q)
  # populate normals
  for(ii in 1:q) {
    XL[ii,ii] <- sqrt(rchisq(1, df = nu-ii+1))
    if(ii > 1) {
      for(jj in 1:(ii-1)) {
        XL[ii,jj] <- rnorm(1)
      }
    }
  }
  if(!inverse) {
    PsiL <- t(chol(Psi))
    XL <- PsiL %*% XL
    X <- XL %*% t(XL)
  } else {
    PsiL <- solve(t(chol(solveV(Psi))))
    XL <- solve(PsiL, XL)
    X <- solve(XL %*% t(XL))
  }
  X
}

# density of matrix normal
dMNormR <- function(X, Lambda, SigmaU, SigmaV, log = FALSE) {
  X <- as.matrix(X)
  Lambda <- as.matrix(Lambda)
  SigmaU <- as.matrix(SigmaU)
  SigmaV <- as.matrix(SigmaV)
  n <- nrow(X)
  p <- ncol(X)
  ans <- n*p*log(2*pi) + sum(diag(solveV(SigmaV, t(X-Lambda)) %*% solveV(SigmaU, X-Lambda)))
  ans <- ans + n*ldetV(SigmaV) + p*ldetV(SigmaU)
  ans <- -ans/2
  if(!log) ans <- exp(ans)
  ans
}

# simulation of matrix normal
rMNormR <- function(Lambda, SigmaU, SigmaV) {
  p <- nrow(Lambda)
  q <- ncol(Lambda)
  Z <- t(rMnorm(q, p))
  X <- Z %*% chol(SigmaV)
  X <- t(chol(SigmaU)) %*% X
  X + Lambda
}

# density of matrix t
dMTR <- function(X, Lambda, SigmaU, SigmaV, nu, log = FALSE) {
  X <- as.matrix(X)
  Lambda <- as.matrix(Lambda)
  SigmaU <- as.matrix(SigmaU)
  SigmaV <- as.matrix(SigmaV)
  n <- nrow(X)
  p <- ncol(X)
  Z <- X - Lambda
  xi <- nu+n+p-1
  ans <- diag(n) + solveV(SigmaU, Z) %*% solveV(SigmaV, t(Z))
  ## .5 * ldet(ans)
  ans <- xi * ldet(ans) + n * ldetV(SigmaV) + p * ldetV(SigmaU) + n*p * log(pi)
  ans <- -.5 * ans + lmgamma(.5*xi, p) - lmgamma(.5*(xi-n), p)
  ## ans <- ans + n*p * log(pi) + lmgamma(.5*(xi-n), p) - lmgamma(.5*xi, p)
  ## ans <- -2 * ans
  if(!log) ans <- exp(ans)
  ans
}

# simulation of matrix T
rMTR <- function(Lambda, SigmaU, SigmaV, nu, prec = FALSE) {
  p <- nrow(Lambda)
  q <- ncol(Lambda)
  V <- rwishR(Psi = SigmaV, nu = nu, inverse = TRUE)
  CL <- t(chol(solveV(V)))
  Z <- t(rMnorm(q, p))
  Z <- Z %*% solve(CL)
  if(!prec) {
    X <- t(chol(SigmaU)) %*% Z + Lambda
  } else {
    X <- solve(chol(SigmaU)) %*% Z + Lambda
  }
  X
}

# density of mniw distribution
dmniwR <- function(X, V, Lambda, Sigma, Psi, nu, log = FALSE) {
  ans <- diwish(V = V, Psi = Psi, nu = nu, log = TRUE)
  if(!all(Sigma == 0)) ans <- ans + dMnorm(X, Lambda, Sigma, V, log = TRUE)
  if(!log) ans <- exp(ans)
  ans
}

# simulation of mniw.
# row variance can be on precision scale by seting prec = TRUE.
rmniwR <- function(Lambda, Sigma, Psi, nu, prec = FALSE) {
  p <- nrow(Lambda)
  q <- ncol(Lambda)
  V <- rwishR(Psi = Psi, nu = nu, inverse = TRUE)
  CL <- t(chol(solveV(V)))
  Z <- t(rMnorm(q, p))
  Z <- Z %*% solve(CL)
  if(!prec) {
    X <- t(chol(Sigma)) %*% Z + Lambda
  } else {
    X <- solve(chol(Sigma)) %*% Z + Lambda
  }
  list(X = X, V = V)
  #PsiL <- solve(t(chol(solveV(Psi))))
  #q <- nrow(PsiL)
  #for(ii in 1:q) {
  #  XL[ii,ii] <- sqrt(rchisq(1, df = nu-ii+1))
  #  if(ii > 1) {
  #    for(jj in 1:(ii-1)) {
  #      XL[ii,jj] <- rnorm(1)
  #    }
  #  }
  #}
  #XL <- solve(XiL, XL)
}

# random effects normal simulation in R
rmNormRER <- function(y, V, lambda, A) {
  C <- solveV(V)
  Q <- solveV(A)
  G <- C + Q
  GU <- chol(G)
  z <- rnorm(length(y))
  mu <- backsolve(r = GU, x = Q %*% (lambda - y), transpose = TRUE)
  drop(backsolve(r = GU, x = z + mu) + y)
}

# multivariate normal simulation in R
rmNormR <- function(mu, V) {
  z <- rnorm(length(mu))
  c(t(chol(V)) %*% z) + mu
}

# multivariate normal density in R
dmNormR <- function(x, mu, V, log = FALSE) {
  p <- length(x)
  ans <- p*log(2*pi) + ldetV(V)
  ans <- ans + crossprod((x-mu), solveV(V, x-mu))[1]
  ans <- -ans/2
  if(!log) ans <- exp(ans)
  ans
}

# unlist to matrices/vectors/scalars to array/matrix/vector.
# dropLast indicates that the last dimension = 1 should be dropped.
unlistMV <- function(X, vtype, dropLast) {
  vtype <- match.arg(vtype, c("matrix", "vector", "scalar"))
  if(vtype == "matrix") {
    PQ <- dim(X[[1]])
    X <- array(unlist(X), dim = c(PQ, length(X)))
    if(dropLast && dim(X)[3] == 1) X <- array(X, dim = PQ)
  } else if(vtype == "vector") {
    X <- do.call(rbind, X)
  } else X <- unlist(X)
  X
}

# create n random p x q matrices.
# If only one of p or q provided the random matrix is a variance.
# rtype is one of "none", "single", "multi".
# if "none", defaults to a matrix of zeros or identity.
# if "single", keep only first entry
# vtype is one of "matrix", "vector", "scalar"
rMM <- function(n, p, q, rtype, vtype) {
  rtype <- match.arg(rtype, c("none", "single", "multi"))
  vtype <- match.arg(vtype, c("matrix", "vector", "scalar"))
  varX <- missing(p) || missing(q) # is X a variance matrix?
  if(rtype != "multi") n <- 1
  if(missing(p)) p <- q
  if(missing(q)) q <- p
  if(vtype == "scalar") {
    X <- lapply(runif(n, 2*q, 3*q),c)
  } else {
    if(rtype == "none") {
      if(!varX) {
        X <- replicate(n, matrix(0,p,q), simplify = FALSE)
      } else {
        X <- replicate(n, diag(p), simplify = FALSE)
      }
    } else {
      if(!varX) {
        X <- replicate(n, rMnorm(p,q), simplify = FALSE)
      } else {
        X <- replicate(n, crossprod(rMnorm(p)) + 5 * diag(p), simplify = FALSE)
      }
    }
    if(vtype == "vector") X <- lapply(X, drop)
  }
  ## if(type != "multi") X <- X[1]
  X
}

# lower of maximum abs/rel error
arDiff <- function(x0, xhat) {
  mx <- abs(xhat - x0)
  min(max(mx), max(mx/abs(x0)))
}

# solve for variance matrices
solveV <- function(V, x, ldV = FALSE) {
  C <- chol(V)
  if (missing(x)) x <- diag(nrow(V))
  ans <- backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
  if (ldV) {
    ldV <- 2 * sum(log(diag(C)))
    ans <- list(y = ans, ldV = ldV)
  }
  ans
}

# check R vs cpp code
expect_Rcpp_equal <- function(fun, icase, mx, ...) {
  Rcpp_comp <- function(fun, icase) mx
  eval(bquote({
    expect_equal(object = Rcpp_comp(.(fun), .(icase)),
                 expected = .(rep(0, length(mx))), ...)
  }))
}

# automatically deduce R-style and cpp-style function arguments from
# case.par specification
# p,q, and drop are common to all the remaining arguments,
# the names of which can be deduced.
# R arguments consist of lists of single arguments, repeated if necessary.
# cpp arguments have appropriate dimensions, and are dropped if necessary.
# args is a named list with named elements p, q, and
# type %in% c("none", "single", "multi")
# drop same as in case.par
get_args <- function(n, args, drop) {
  anames <- names(args)
  # simulate argument values
  sargs <- sapply(args, function(vargs) {
    do.call(rMM, c(list(n = n), vargs))
  }, simplify = FALSE)
  # R format
  rargs <- sapply(anames, function(nm) {
    if(args[[nm]]$rtype != "multi") {
      rep(sargs[[nm]], n)
    } else sargs[[nm]]
  }, simplify = FALSE)
  # cpp format
  cargs <- sapply(anames, function(nm) {
    unlistMV(sargs[[nm]],
             vtype = args[[nm]]$vtype, drop = drop)
  }, simplify = FALSE)
  cargs <- cargs[sapply(args, function(aa) aa$rtype != "none")]
  list(R = rargs, cpp = cargs)
}

# determine if all arguments are single
all_single <- function(cp) {
  nignore <- c("drop", "inverse", "prec", "p", "q", "r")
  single <- sapply(cp[!names(cp) %in% nignore],
                   function(x) {
                     if(!x %in% c("none", "single", "multi")) {
                       stop("wrong input to all_single.")
                     }
                     x != "multi"
                   })
  all(single)
}
