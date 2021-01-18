\name{pmvn}
\alias{pmvn}
\title{Quasi-Monte Carlo method for multivariate normal probabilities}
\description{
  Compute multivariate normal probabilities with the dense-matrix based Quasi-Monte Carlo method and the tile-low-rank-matrix based Quasi-Monte Carlo method.
}
\usage{
pmvn(lower = -Inf, upper = Inf, mean = 0, sigma = NULL, 
uselog2 = FALSE, algorithm = GenzBretz(), ...)
}
\arguments{
  \item{lower}{lower integration limits, a numeric vector of length n}
  \item{upper}{upper integration limits, a numeric vector of length n}
  \item{mean}{the mean parameter, a numeric vector of length n}
  \item{sigma}{the covariance matrix of dimension n}
  \item{uselog2}{whether return the result as the logarithm to the base 2}
  \item{algorithm}{an object of class \code{\link{GenzBretz}} or
    \code{\link{TLRQMC}} defining the hyper parameters of this algorithm}
  \item{...}{additional parameters used to construct `sigma` when it is not given:
    \itemize{
      \item{geom}{a matrix of dimension n-by-2, specifying n spatial locations in the 2D unit square}
      \item{type}{the name of the covariance kernel, a string. Currently, only the Matern covariance function, e.g., "matern", is supported. Not case-sensitive. It should be given when `sigma` is not given}
      \item{para}{the parameter for the covariance kernel, a numeric vector. When `type` is "matern", the length of `para` should be 4, representing the scale, range, smoothness, and nugget parameters of the covariance function. It should be given when `sigma` is not given}
    }
  }
}
\value{
  When `uselog2` is set FALSE, the function returns the estimated probability with one attribute of the estimation error. 
  When `uselog2` is set TRUE, the function only returns the estimated log-probability to the base 2. This is useful when the estimated probability is smaller than the machine precision.
}
\details{
  When `algorithm` is of the class `GenzBretz`, the Quasi-Monte Carlo sampling described in Genz, A. (1992) is used.
  When `algorithm` is of the class `TLRQMC`, the Quasi-Monte Carlo sampling with the tile-low-rank representation of the covariance matrix, described in Cao et al. (2020), is used.
  When `sigma`, is given, `geom`, `type`, and `para` are not used. Otherwise, a covariance matrix is created with the information from `geom`, `type`, and `para`.
}
\references{
  Genz, A. (1992), "Numerical computation of multivariate
normal probabilities," Journal of Computational and
Graphical Statistics, 1, 141-149.
  Cao, J., Genton, M. G., Keyes, D. E., & Turkiyyah, G. M. (2019). "Exploiting Low
Rank Covariance Structures for Computing High-Dimensional Normal and Student-
t Probabilities"
}
\examples{
  n = 225
  set.seed(0)
  a = rep(-10, n)
  b = rnorm(n, 3, 2)
  m = 15
  epsl = 1e-4
  vec1 = 1 : m
  vec2 = rep(1, m)
  geom = cbind(kronecker(vec1, vec2), kronecker(vec2, vec1))
  geom = geom / m
  beta = 0.3
  idx = zorder(geom)
  geom = geom[idx, ]
  a = a[idx]
  b = b[idx]
  distM = as.matrix(dist(geom))
  covM = exp(-distM / beta)
  pmvn(lower = a, upper = b, mean = 2, sigma = covM, uselog2 = FALSE, 
       algorithm = GenzBretz(N = 521))
  pmvn(lower = a, upper = b, mean = 2, uselog2 = TRUE, geom = geom, 
       kernelType = "matern", para = c(1.0, 0.3, 0.5, 0.0))
  pmvn(lower = a, upper = b, mean = 2, sigma = covM, uselog2 = FALSE, 
       algorithm = TLRQMC(N = 521, m = m, epsl = epsl))
  pmvn(lower = a, upper = b, mean = 2, uselog2 = TRUE, geom = geom, 
       algorithm = TLRQMC(N = 521, m = m, epsl = epsl),
       kernelType = "matern", para = c(1.0, 0.3, 0.5, 0.0))
}
\author{Jian Cao, Marc Genton, David Keyes, George Turkiyyah}

