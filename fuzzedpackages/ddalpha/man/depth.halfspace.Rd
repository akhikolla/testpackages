\name{depth.halfspace}
\alias{depth.halfspace}
\title{
Calculate Halfspace Depth
}
\description{
Calculates the exact or random Tukey (=halfspace, location) depth  (Tukey, 1975) of points w.r.t. a multivariate data set.
}
\usage{
depth.halfspace(x, data, exact, method, num.directions = 1000, seed = 0)
}
\arguments{
  \item{x}{
Matrix of objects (numerical vector as one object) whose depth is to be calculated; each row contains a \eqn{d}-variate point. Should have the same dimension as \code{data}.
}
  \item{data}{
Matrix of data where each row contains a \eqn{d}-variate point, w.r.t. which the depth is to be calculated.
}
  \item{exact}{
The type of the used method. The default is \code{exact=F}, which leads to approximate computation of the Tukey depth. For \code{exact=F}, \code{method="Sunif.1D"} is used by default. If \code{exact=T}, the Tukey depth is computed exactly, with \code{method="recursive"} by default.}
  \item{method}{
For \code{exact=F}, if \code{method="Sunif.1D"} (by default), the Tukey depth is computed approximately by being minimized over univariate projections (see Details below).

For \code{exact=T}, the Tukey depth is calculated as the minimum over all combinations of \eqn{k} points from \code{data} (see Details below). In this case parameter \code{method} specifies \eqn{k}, with possible values \eqn{1} for \code{method="recursive"} (by default), \eqn{d-2} for \code{method="plane"}, \eqn{d-1} for \code{method="line"}.

The name of the method may be given as well as just parameter \code{exact}, in which case the default method will be used.
}
  \item{num.directions}{
Number of random directions to be generated (for \code{method="Sunif.1D"}). The algorithmic complexity is linear in the number of observations in \code{data}, given the number of directions.
}
  \item{seed}{
The random seed. The default value \code{seed=0} makes no changes (for \code{method="Sunif.1D"}).
}
}
\details{
For \code{exact=F}, if \code{method="Sunif.1D"}, the Tukey depth is computed approximately using the random Tukey depth method proposed by Cuesta-Albertos and Nieto-Reyes (2008). Here the depth is determined as the minimum univariate Tukey depth of the - on lines in several directions - projected data. The directions are distributed uniformly on the \eqn{(d-1)}-sphere; the same direction set is used for all points.

For \code{exact=T}, the Tukey depth is computed exactly as the minimum of the sum of the depths in two orthogonal complementary affine subspaces, which dimensions add to \eqn{d}: one of the subspaces (combinatorial) is the \eqn{k}-dimensional hyperplane through (a point from) \code{x} and \eqn{k} points from \code{data}, another one is its orthogonal complement (see Dyckerhoff and Mozharovskyi, 2016 for the detailed description of the algorithmic framework). The algorithm then minimizes the depth over all combinations of \eqn{k} points, in which the depth in the orthogonal complements is computed using an exact algorithm. In this case, parameter \code{method} specifies the dimensionality \eqn{k} of the combinatorial space. The implemented (reasonable) algorithms (and corresponding names) are: \eqn{k=1} (or \code{method="recursive"}), \eqn{k=d-2} (or \code{method="plane"}), and \eqn{k=d-1} (or \code{method="line"}).
}
\value{
Numerical vector of depths, one for each row in \code{x}; or one depth value if \code{x} is a numerical vector.
}
\references{
Cuesta-Albertos, J.A. and Nieto-Reyes, A. (2008). The random Tukey depth. \emph{Computational Statistics and Data Analysis} \bold{52} 4979--4988.

Dyckerhoff, R. and Mozharovskyi, P. (2016). Exact computation of the halfspace depth. \emph{Computational Statistics and Data Analysis} \bold{98} 19--30.

Rousseeuw, P.J. and Ruts, I. (1996). Algorithm AS 307: Bivariate location depth. \emph{Journal of the Royal Statistical Society. Seriec C (Applied Statistics)} \bold{45} 516--526.

Tukey, J.W. (1974). Mathematics and the picturing of data. In: \emph{Proceeding of the International Congress of Mathematicians}, Vancouver, 523--531.
}
\seealso{
\code{\link{depth.Mahalanobis}} for calculation of Mahalanobis depth.

\code{\link{depth.projection}} for calculation of projection depth.

\code{\link{depth.simplicial}} for calculation of simplicial depth.

\code{\link{depth.simplicialVolume}} for calculation of simplicial volume depth.

\code{\link{depth.spatial}} for calculation of spatial depth.

\code{\link{depth.zonoid}} for calculation of zonoid depth.

\code{\link{depth.potential}} for calculation of data potential.

}
\examples{
# 3-dimensional normal distribution
data <- mvrnorm(200, rep(0, 3), 
                matrix(c(1, 0, 0,
                         0, 2, 0, 
                         0, 0, 1),
                nrow = 3))
x <- mvrnorm(10, rep(1, 3), 
             matrix(c(1, 0, 0,
                      0, 1, 0, 
                      0, 0, 1),
             nrow = 3))
              
# default - random Tukey depth
depths <- depth.halfspace(x, data)
cat("Depths: ", depths, "\n")

# default exact method - "recursive"
depths <- depth.halfspace(x, data, exact = TRUE)
cat("Depths: ", depths, "\n")

# method "line"
depths <- depth.halfspace(x, data, method = "line")
cat("Depths: ", depths, "\n")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
