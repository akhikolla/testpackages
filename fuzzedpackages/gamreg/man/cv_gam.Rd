\name{cv.gam}
\alias{cv.gam}
\title{Robust Cross-Validation}
\usage{
  cv.gam(X, Y, init.mode = c("sLTS", "RLARS", "RANSAC"),
         lambda.mode = "lambda0", lmax = 1, lmin = 0.05, nlambda = 50,
         fold = 10, ncores = 1, gam = 0.1, gam0 = 0.5, intercept = "TRUE",
         alpha = 1, ini.subsamp = 0.2, ini.cand = 1000, alpha.LTS = 0.75,
         nlambda.LTS = 40)
}
\description{
  Compute Robust Cross-Validation for selecting best model.
}
\arguments{
  \item{X}{Predictor variables Matrix.}
  \item{Y}{Response variables Matrix.}
  \item{init.mode}{\code{"sLTS"}: a initial point is the estimate of sparse least trimmed squares. \code{"RLARS"}: a initial point is the estimate of Robust LARS. \code{"RANSAC"}: a initial point is the estimate of RANSAC algorithm.}
  \item{lambda.mode}{\code{"lambda0"}: Robust Cross-Validation uses grids on range [0.05lambda0,lambda0] with log scale, where lambda0 is an estimator of sparse tuning parameter which would shrink regression coefficients to zero.}
  \item{lmax}{When \code{lambda.mode} is not lambda0, upper bound of range of grids is lmax.}
  \item{lmin}{When \code{lambda.mode} is not lambda0, lower bound of range of grids is lmin.}
  \item{nlambda}{The number of grids for Robust Cross-Validation.}
  \item{fold}{the number of folds for K-fold Robust Cross-Validation. If \code{fold} equals to sample size, Robust Cross-Validation is leave-one-out method.}
  \item{ncores}{positive integer giving the number of processor cores to be used for parallel computing (the default is 1 for no parallelization).}
  \item{gam}{Robust tuning parameter of gamma-divergence for regression.}
  \item{gam0}{tuning parameter of Robust Cross-Validation.}
  \item{intercept}{Should intercept be fitted \code{TRUE} or set to zero \code{FALSE}}
  \item{alpha}{The elasticnet mixing parameter, with \eqn{0 \le \alpha \le 1}. \code{alpha=1} is the lasso penalty, and \code{alpha=0} the ridge penalty.}
  \item{ini.subsamp}{The fraction of subsamples in "\code{RANSAC}".}
  \item{ini.cand}{The number of candidates for estimating itnial points in "\code{RANSAC}".}
  \item{alpha.LTS}{The fraction of subsamples for trimmed squares in "\code{sLTS}".}
  \item{nlambda.LTS}{The number of grids for sparse tuning parameter in "\code{sLTS}".}
}
\value{
  \item{lambda}{A numeric vector giving the values of the penalty parameter.}
  \item{fit}{All results at each lambda.}
  \item{Rocv}{The result of best model by Robust Cross-Validation.}
}
\references{
  Kawashima, T. and Fujisawa, H. (2017).
  \emph{Robust and Sparse Regression via gamma-divergence, Entropy, 19(11).}\cr
  Fujisawa, H. and Eguchi, S. (2008).
  \emph{Robust parameter estimation with a small bias against heavy contamination, Journal of Multivariate Analysis, 99(9), 2053-2081.}\cr
}

\details{
  If the "\code{RANSAC}" is used as the initial point, the parameter \code{ini.subsamp} and \code{ini.cand} can be determined carefully. The smaller \code{ini.subsamp} is, the more robust initial point is. However, less efficiency.
}
\examples{
  ## generate data
  library(mvtnorm)
  n <- 30                      # number of observations
  p <- 10                      # number of expalanatory variables

  epsilon <- 0.1               # contamination ratio

  beta0 <- 0.0                 # intercept
  beta <- c(numeric(p))        # regression coefficients
  beta[1] <- 1
  beta[2] <- 2
  beta[3] <- 3
  beta[4] <- 4

  Sigma <- 0.2^t(sapply(1:p, function(i, j) abs(i-j), 1:p))
  X <- rmvnorm(n, sigma=Sigma) # explanatory variables
  e <- rnorm(n) # error terms

  i <- 1:ceiling(epsilon*n)    # index of outliers
  e[i] <- e[i] + 20            # vertical outliers
  Y <- beta0*(numeric(n)+1) + X\%*\%beta


  res <- cv.gam(X,Y,nlambda = 5, nlambda.LTS=20 ,init.mode="sLTS")


}
\author{
  Takayuki Kawashima
}
