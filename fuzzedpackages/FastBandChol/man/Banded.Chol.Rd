
\name{banded.chol}
\alias{banded.chol}
\title{Computes estimate of covariance matrix by banding the Cholesky factor}
\usage{
banded.chol(X, bandwidth, centered = FALSE)
}
\arguments{
\item{X}{A data matrix with \eqn{n} rows and \eqn{p} columns. Rows are assumed to be independent realizations from a \eqn{p}-variate distribution with covariance \eqn{\Sigma}.}

\item{bandwidth}{A positive integer. Must be less than \eqn{n-1} and \eqn{p-1}.}

\item{centered}{Logical. Is data matrix centered? Default is \code{centered = FALSE}}
}
\value{
A list with
\item{est}{The estimated covariance matrix.}
}
\description{
Computes estimate of covariance matrix by banding the Cholesky factor using a modified Gram Schmidt algorithm implemented in RcppArmadillo. 
}
\examples{
## set sample size and dimension
n=20
p=100

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

## compute estimate
out1 = banded.chol(X, bandwidth=4)
}

