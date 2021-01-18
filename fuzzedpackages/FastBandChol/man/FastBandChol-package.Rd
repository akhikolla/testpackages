\name{FastBandChol-package}
\alias{FastBandChol-package}
\alias{FastBandChol}
\docType{package}
\title{
Fast estimation of covariance matrix by banded Cholesky factor
}
\description{
Fast and numerically stable estimation of covariance matrix by banding the Cholesky factor using a modified Gram-Schmidt algorithm implemented in RcppArmadilo. See <https://stat.umn.edu/~molst029> for details on the algorithm. 
}
\details{
\tabular{ll}{
Package: \tab FastBandChol\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2015-08-22\cr
License: \tab GPL-2\cr
}

}
\author{
Aaron Molstad
}
\references{
Rothman, A.J., Levina, E., and Zhu, J. (2010). A new approach to Cholesky-based covariance regularization in high dimensions. Biometrika, 97(3):539-550.
}

\examples{
## set sample size and dimension
n = 20
p = 100

## create covariance with AR1 structure
Sigma = matrix(0, nrow=p, ncol=p)
for(l in 1:p){
  for(m in 1:p){
    Sigma[l,m] = .5^(abs(l-m))
  }
}

## simulation Normal data
eo1 = eigen(Sigma)
Sigma.sqrt = eo1$vec\%*\%diag(eo1$val^.5)\%*\%t(eo1$vec)
X = t(Sigma.sqrt\%*\%matrix(rnorm(n*p), nrow=p, ncol=n))

## compute estimates
est.sample = banded.sample(X, bandwidth=4)$est
est.chol = banded.chol(X, bandwidth=4)$est


}
