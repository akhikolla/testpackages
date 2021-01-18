\name{unilps}
\alias{unilps}
\title{
  univariate model fitting with lasso penalty
}
\description{
  fit univariate log-linear model using accelerated block-coordinate
  relaxation algorithm.
}
\usage{
unilps(formula, data = list(),
       family = c("gaussian", "binomial"),
       lambda = NULL, nlambda = 100,
       lambda.min.ratio = ifelse(nobs<nvars, .01, .0001),
       output = 0, tune = c("AIC", "BIC", "GACV", "BGACV"))
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{an optional data frame containing the variables in the
    model. By default the variables are taken from the
    environment from which \code{unifit} is called.}
  \item{family}{
    a GLM family, currently suport gaussian and binomial.
  }
  \item{lambda}{
    a user specified tuning sequece. Typical usage is to have the
    program compute its own \code{lambda}.}
  \item{nlambda}{the number of \code{lambda} values, default is 100.}
  \item{lambda.min.ratio}{Smallest value for \code{lambda}, as a
    fraction of \code{lambda.max}. The default depends on the sample
    size \code{nobs} relative to the number of variables.}
  \item{output}{with values 0 or 1, indicating whether the fitting
    process is muted or not.}
  \item{tune}{tuning approach, available methods including AIC, BIC,
    GACV, BGACV.}
}
\details{
  The \code{unilps} utilize the class structure of the underlying C++
  code and fitted the model with accelerated block-coordinate relaxation algorithm.
}
\value{
  An object of classes \code{mvbfit} and \code{lps}, for which some methods are
  available.
}
\seealso{
  \code{unilps}, \code{mvblps}
}
\examples{
n <- 100
p <- 4
x <- matrix(rnorm(n * p, 0, 4), n, p)
eta <- x %*% c(2, -1, rep(0, p - 2))
pr <- exp(eta) / (1+ exp(eta))
res <- rbinom(n, 1, pr)
fit <- unilps(res ~ x - 1, family = 'binomial')
}