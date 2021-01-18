\name{loglike}
\alias{loglike}
\title{
  negative loglikelihood evaluation
}
\description{
  evaluate negative loglikelihood of the corresponding family of model.
}
\usage{
loglike(x, y, input,
        family = c("gaussian", "bernoulli", "mvbernoulli")) 
}
\arguments{
  \item{x}{design matrix.}
  \item{y}{output binary matrix with number of columns equal to the
    number of outcomes per observation.}
  \item{input}{vector of the fitted coefficients for the distribution family.}
  \item{family}{
    a GLM family, currently support gaussian, binomial and mvbernoulli
    (multivariate Bernoulli).
  }
}
\details{
  evaluate the negative log-likelihood to examine the performance of the model.
}
\value{
  a double value returned as the negative log-likelihood
}
\seealso{
  \code{unifit}, \code{mvbfit}
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
loglike(x, res$response, fitMVB$beta, "mvbernoulli")
}