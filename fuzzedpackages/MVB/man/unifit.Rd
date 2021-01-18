\name{unifit}
\alias{unifit}
\title{
  univariate model fitting
}
\description{
  fit univariate log-linear model using Newton-Raphson algorithm.
}
\usage{
unifit(formula, data = list(),
       family = c("gaussian", "binomial"),
       output = 0)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{an optional data frame containing the variables in the
    model. By default the variables are taken from the
    environment from which \code{unifit} is called.}
  \item{family}{
    a GLM family, currently support gaussian and binomial.
  }
  \item{output}{with values 0 or 1, indicating whether the fitting
    process is muted or not.}
}
\details{
  The \code{unifit} utilize the class structure of the underlying C++
  code and fitted the model with Newton-Raphson algorithm.
}
\value{
  An object of class \code{mvbfit}, for which some methods are
  available.
}
\seealso{
  \code{unilps}, \code{mvbfit}
}
\examples{
n <- 100
p <- 4
x <- matrix(rnorm(n * p, 0, 4), n, p)
eta <- x %*% c(2, -1, rep(0, p - 2))
pr <- exp(eta) / (1+ exp(eta))
res <- rbinom(n, 1, pr)
fit <- unifit(res ~ x - 1, family = 'binomial')
}