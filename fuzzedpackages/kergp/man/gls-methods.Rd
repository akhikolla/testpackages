\name{gls-methods}
\docType{methods}
\alias{gls,covAll-method}
\alias{gls-methods}
\title{
  Generalized Least Squares Estimation with a Given Covariance Kernel
}

\description{
  
  Generalized Least Squares (GLS) estimation for a linear model with
  a covariance given by the covariance kernel object. The method gives
  auxiliary variables as needed in many algebraic computations.

}
\details{

  There are two options: for unknown trend, this is the usual GLS
  estimation with given covariance kernel; for a known trend, it returns
  the corresponding auxiliary variables (see \code{value} below).

}

\usage{

\S4method{gls}{covAll}(object,
    y, X, F = NULL, varNoise = NULL, 
    beta = NULL, checkNames = TRUE,
    ...)
}

\arguments{
  \item{object}{

    An object with \code{"covAll"} class.

  }
  \item{y}{

    The response vector with length \eqn{n}.

  }
  \item{X}{

    The input (or spatial design) matrix with \eqn{n} rows and \eqn{d}
    columns. This matrix must be compatible with the given covariance
    object, see \code{\link{checkX,covAll,matrix-method}}.

  }
  \item{F}{

    A trend design matrix with \eqn{n} rows and \eqn{p} columns. When
    \code{F} is \code{NULL} no trend is used and the response \code{y}
    is simply a realization of a centered Gaussian Process with covariance kernel given by \code{object}.

  }
  \item{varNoise}{

    A known noise variance. When provided, must be a positive numeric
    value.

  }
  \item{beta}{

    A known vector of trend parameters. Default is \code{NULL}
    indicating that the trend parameters must be estimated.

  }
  \item{checkNames}{

    Logical. If \code{TRUE} (default), check the compatibility of
    \code{X} with \code{object}, see \code{\link{checkX}}.

  }
  \item{\dots}{

    not used yet.

  }
}

% \details{
%  This function computes the estimate of a regression with
%  covariance structure given by \code{object} with possibly extra
%  noise.
%
%  The regression is transformed to an OLS regression using left
%  multiplication by the inverse (Lower) Cholesky root
%  \eqn{\mathbf{L}^{-1}}{inv(L)} of the covariance matrix
%  \eqn{\mathbf{C}}{C}. 
%  Then a classical QR decomposition of the (transformed) design matrix
%  is used to get the estimated parameter and error.
%}

\value{
  
  A list with several elements.
  \item{betaHat}{

    Vector \eqn{\widehat{\boldsymbol{\beta}}}{} of length \eqn{p}
    containing the estimated coefficients if \code{beta = NULL}, or the
    known coefficients \eqn{\boldsymbol{\beta}}{beta} either.

  }
  \item{L}{

    The (lower) Cholesky root matrix \eqn{\mathbf{L}}{L} of the
    covariance matrix \eqn{\mathbf{C}}{C}. This
    matrix has \eqn{n} rows and \eqn{n} columns and \eqn{\mathbf{C} =
    \mathbf{L} \mathbf{L}^\top}{C = L \%*\% t(L)}.

  }
  \item{eStar}{

    Vector of length \eqn{n}: \eqn{\mathbf{e}^\star = \mathbf{L}^{-1}(\mathbf{y} -
    \mathbf{X} \widehat{\boldsymbol{\beta}})}{eStar = inv(L) \%*\% (y - X betaHat)}.

  }
  \item{Fstar}{

    Matrix \eqn{n \times p}{(n, p)}:  \eqn{\mathbf{F}^\star :=
    \mathbf{L}^{-1}\mathbf{F}}{FStar = inv(L) \%*\% F}.

  } 
  \item{sseStar}{

    Sum of squared errors:
    \eqn{{\mathbf{e}^\star}^\top\mathbf{e}^\star}{t(eStar) \%*\% eStar}.

  }
  \item{RStar}{

    The 'R' upper triangular \eqn{p \times p}{(p, p)} matrix in the QR
    decomposition of \code{FStar}: \eqn{\mathbf{F}^\star =
      \mathbf{Q}\mathbf{R}^\star}{FStar = Q \%*\% RStar}.
      

  }

  All objects having length \eqn{p} or having one of their dimension
  equal to \eqn{p} will be \code{NULL} when \code{F} is \code{NULL},
  meaning that \eqn{p = 0}.
  
}
\references{

  Kenneth Lange (2010), \emph{Numerical Analysis for Statisticians} 2nd ed.
  pp. 102-103. Springer-Verlag,

}

\author{
  
  Y. Deville, O. Roustant

}

\examples{
## a possible 'covTS'
myCov <- covTS(inputs = c("Temp", "Humid"),
               kernel = "k1Matern5_2",
               dep = c(range = "input"),
               value = c(range = 0.4))
d <- myCov@d; n <- 100;
X <- matrix(runif(n*d), nrow = n, ncol = d)
colnames(X) <- inputNames(myCov)
## generate the 'GP part'  
C <- covMat(myCov, X = X)
L <- t(chol(C))
zeta <- L \%*\% rnorm(n)
## trend matrix 'F' for Ordinary Kriging
F <- matrix(1, nrow = n, ncol = 1)
varNoise <- 0.5
epsilon <- rnorm(n, sd = sqrt(varNoise))
beta <- 10
y <- F \%*\% beta + zeta + epsilon
fit <- gls(myCov, X = X, y = y, F = F, varNoise = varNoise)

}
