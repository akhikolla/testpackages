\name{stepfit}
\alias{stepfit}
\title{
  step-wisd multivariate model fitting
}
\description{
  stepwise fit multivariate log-linear Bernoulli model using Newton-Raphson algorithm.
}
\usage{
stepfit(x, y, maxOrder = 2,
        output = 0,
        direction = c("backward", "forward"),
        tune = c("AIC", "BIC", "GACV", "BGACV"),
        start = NULL)
}
\arguments{
  \item{x}{input design matrix.}
  \item{y}{output binary matrix with number of columns equal to the
    number of outcomes per observation.}
  \item{maxOrder}{maximum order of interactions to be considered in outcomes.}
  \item{output}{with values 0 or 1, indicating whether the fitting
    process is muted or not.}
  \item{direction}{the mode of stepwise search and default is backward.}
  \item{tune}{tuning approach, available methods including AIC, BIC,
    GACV, BGACV.}
  \item{start}{starting object of type mvbfit.}
}
\details{
  The \code{stepfit} utilize the class structure of the underlying C++
  code and stepwisd fitted the model with Newton-Raphson algorithm.
}
\value{
  An object of class \code{mvbfit}, for which some methods are
  available.
}
\seealso{
  \code{mvblps}, \code{unifit}, \code{stepfit}, \code{mvb.simu}
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