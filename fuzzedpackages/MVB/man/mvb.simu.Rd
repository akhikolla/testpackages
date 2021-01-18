\name{mvb.simu}
\alias{mvb.simu}
\title{
  generate multivariate Bernoulli simulated data
}
\description{
  for given coefficients and design matrix, generate the corresponding
  responses according multivariate Bernoulli model 
}
\usage{
mvb.simu(coefficients, x, K = 2, offset = as.double(0))
}
\arguments{
  \item{coefficients}{coefficients matrix, number of columns should be
    less than \code{2^K}.}
  \item{x}{design matrix.}
  \item{K}{number of outcomes for the model.}
  \item{offset}{non-penalized terms in coefficients, corresponding to a
    unit column in design matrix, which is generated automaticly.}
}
\details{
  The response variables are simulated according to cononical link
  function of multivariate Bernoulli model with coefficients speicified.
}
\value{
  \item{response}{matrix for outcomes, with dimension \code{nobs} times
    \code{K}.}
  \item{beta}{expanded coefficients from input argument
    \code{coefficients} and \code{offset}.}
}
\seealso{
  \code{mvbfit}, \code{mvblps}
}
\examples{
# fit a simple MVB log-linear model
n <- 1000
p <- 5
kk <- 2
tt <- NULL
alter <- 1
for (i in 1:kk) {
  vec <- rep(0, p)
  vec[i] <- alter
  alter <- alter * (-1)
  tt <- cbind(tt, vec)
}
tt <- 1.5 * tt
tt <- cbind(tt, c(rep(0, p - 1), 1))

x <- matrix(rnorm(n * p, 0, 4), n, p)
res <- mvb.simu(tt, x, K = kk, rep(.5, 2))
fitMVB <- mvbfit(x, res$response, output = 1)
}