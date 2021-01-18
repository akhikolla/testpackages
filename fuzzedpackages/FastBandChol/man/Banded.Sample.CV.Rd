\name{banded.sample.cv}
\alias{banded.sample.cv}
\title{Selects bandwidth for sample covariance matrix by cross validation}
\usage{
banded.sample.cv(X, bandwidth, folds = 3, est.eval = TRUE, Frob = TRUE)
}
\arguments{
\item{X}{A data matrix with \eqn{n} rows and \eqn{p} columns. Rows are assumed to be independent realizations from a \eqn{p}-variate distribution with covariance \eqn{\Sigma}.}

\item{bandwidth}{A vector of candidate bandwidths. Candidate bandwidths can only positive integers such that the maximum is less than \eqn{p-1}}.


\item{folds}{The number of folds used for cross validation. Default is \code{folds =3}.}

\item{est.eval}{Logical: \code{est.eval = TRUE} returns a list with both the selected bandwidth and the estimated covariance matrix. \code{est.eval=FALSE} returns a list with only the selected bandwidth. The default is \code{est.eval = TRUE}.}

\item{Frob}{Logical: \code{Frob = TRUE} uses squared Frobenius norm loss for cross-validation. \code{Frob = FALSE} uses operator norm loss. Default is \code{Frob = TRUE}.}
}

\value{
A list with
\item{bandwidth.min}{the bandwidth minimizing cv error}

\item{est}{the sample covariance matrix at bandwidth.min}
}
\description{
Selects bandwidth for sample covariance matrix by k-fold cross validation
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

## perform cross validation
k = 4:7
out2.cv = banded.sample.cv(X, bandwidth=k, folds=5)
}

