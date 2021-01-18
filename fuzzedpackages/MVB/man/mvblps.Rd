\name{mvblps}
\alias{mvblps}
\title{
  multivariate Bernoulli LASSO model fitting 
}
\description{
  fit multivariate Bernoulli LASSO model accelerated block-coordinate relaxation algorithm.
}
\usage{
mvblps(x, y, maxOrder = 2, lambda = NULL, nlambda = 100,
       lambda.min.ratio = ifelse(nobs<nvars, .01, .0001),
       output = 0, printIter = 100, search = c('nm', 'grid'),
       tune = c("AIC", "BIC", "GACV", "BGACV"))
}
\arguments{
  \item{x}{input design matrix.}
  \item{y}{output binary matrix with number of columns equal to the
    number of outcomes per observation.}
  \item{maxOrder}{maximum order of interactions to be considered in outcomes.}
  \item{lambda}{
    a user specified tuning sequece. Typical usage is to have the
    program compute its own \code{lambda}.}
  \item{nlambda}{the number of \code{lambda} values, default is 100.}
  \item{lambda.min.ratio}{Smallest value for \code{lambda}, as a
    fraction of \code{lambda.max}. The default depends on the sample
    size \code{nobs} relative to the number of variables.}
  \item{output}{with values 0 or 1, indicating whether the fitting
    process is muted or not.}
  \item{printIter}{Number of iterations to be printed if output is
    true.}
  \item{search}{Tuning search approach, \code{nm} for Nelder Mead and
    \code{grid} for grid search.}
  \item{tune}{tuning approach, available methods including AIC, BIC,
    GACV, BGACV.}

}
\details{
  The \code{mvblps} utilize the class structure of the underlying C++
  code and fitted the model with accelerated block-coordinate relaxation algorithm.
}
\value{
  An object of classes \code{mvbfit} and \code{lps}, for which some methods are
  available.
}
\seealso{
  \code{mvbfit}, \code{unifit}, \code{stepfit}, \code{mvb.simu}
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
fitMVB <- mvblps(x, res$response, output = 1)

}