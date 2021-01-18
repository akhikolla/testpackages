\name{ldmvnorm}
\alias{ldmvnorm}
\title{Multivariate Normal Log-Density for Complete and Incomplete Data}
\description{
This function provides the log-density function for the multivariate normal distribution with mean equal to mu and covariance matrix Sigma.
Marginal distributions will be used when the vector (or matrix) of quantiles is incomplete. That is, when the vector (or matrix) of quantiles
contain \code{NA}.
}
\usage{
ldmvnorm(x, mu, Sigma, onNA=0)
}
\arguments{
  \item{x}{Vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.}
  \item{mu}{Mean vector, default is \code{rep(0, length = ncol(x))}}
  \item{Sigma}{Covariance matrix, default is \code{diag(ncol(x))}.}
  \item{onNA}{Action for a row on NAs. Default is to return \code{0}.}
}
\author{Alessio Farcomeni \email{alessio.farcomeni@uniroma1.it}, Andy
Leung \email{andy.leung@stat.ubc.ca}}
\examples{
x <- matrix(rnorm(1000),100, 10)
u <- matrix(rbinom(1000, 1, 0.1), 100, 10)
x[ u == 1 ] <- NA
mu <- rep(0,10)
Sigma <- diag(10)
ldmvnorm(x, mu, Sigma)
}


