\name{pmvt}
\alias{pmvt}
\title{Quasi-Monte Carlo method for Student-$t$ probabilities}
\description{
  Compute multivariate Student-$t$ probabilities with the dense-matrix based Quasi-Monte Carlo method and the tile-low-rank-matrix based Quasi-Monte Carlo method.
}
\usage{
pmvt(lower = -Inf, upper = Inf, delta = 0, df = 1, sigma = NULL, 
     uselog2 = FALSE, algorithm = GenzBretz(), 
     type = "Kshirsagar", ...)
}
\arguments{
  \item{lower}{lower integration limits, a numeric vector of length n}
  \item{upper}{upper integration limits, a numeric vector of length n}
  \item{delta}{the vector of noncentrality parameters of length n, for type = "shifted" delta specifies the mode}
  \item{df}{a positive numeric value denoting the degrees of freedom}
  \item{sigma}{the covariance matrix of dimension n}
  \item{uselog2}{whether return the result as the logarithm to the base 2}
  \item{algorithm}{an object of class \code{\link{GenzBretz}} or \code{\link{TLRQMC}} defining the hyper parameters of this algorithm}
  \item{type}{type of the noncentral multivariate $t$ distribution to be computed. `type` = "Kshirsagar" corresponds to formula (1.4) in Genz and Bretz (2009). `type` = "shifted" corresponds to the formula right before formula (1.4) in Genz and Bretz (2009)}
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
  Genz, A. and Bretz, F. (2009), Computation of Multivariate Normal and t Probabilities. Lecture Notes in Statistics, Vol. 195. Springer-Verlag, Heidelberg.
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
  df = 10
  pmvt(lower = a, upper = b, delta = 2, df = df,
       sigma = covM, uselog2 = FALSE, algorithm = GenzBretz(N = 521),
       type = "Kshirsagar")
  pmvt(lower = a, upper = b, delta = 2, df = df,
       uselog2 = TRUE, type = "shifted", geom = geom, 
       kernelType = "matern", para = c(1.0, 0.3, 0.5, 0.0))
  pmvt(lower = a, upper = b, delta = 2, df = df,
       sigma = covM, uselog2 = FALSE, 
       algorithm = TLRQMC(N = 521, m = m, epsl = epsl), type = "Kshirsagar")
  pmvt(lower = a, upper = b, delta = 2, df = df,
       uselog2 = TRUE, type = "shifted", geom = geom,
       algorithm = TLRQMC(N = 521, m = m, epsl = epsl),
       kernelType = "matern", para = c(1.0, 0.3, 0.5, 0.0))
}
\author{Jian Cao, Marc Genton, David Keyes, George Turkiyyah}

