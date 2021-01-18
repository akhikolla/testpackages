\name{depth.qhpeeling}
\alias{depth.qhpeeling}
\title{
Calculate Convex Hull Peeling Depth
}
\description{
Calculates the convex hull peeling depth of points w.r.t. a multivariate data set.
}
\usage{
depth.qhpeeling(x, data)
}
\arguments{
  \item{x}{
Matrix of objects (numerical vector as one object) whose depth is to be calculated; each row contains a \eqn{d}-variate point. Should have the same dimension as \code{data}.
}
  \item{data}{
Matrix of data where each row contains a \eqn{d}-variate point, w.r.t. which the depth is to be calculated.
}
}
\details{
Calculates the convex hull peeling depth (Eddy, 1982; see also Cascos, 2009).
}
\value{
Numerical vector of depths, one for each row in \code{x}; or one depth value if \code{x} is a numerical vector. Each depth value equals the number of the convex hulls to be peeled from \code{data} so that (the corresponding row of) \code{x} is not contained in the convex hull of the rest of the data; the depths are normalized by the number of points in \code{data}.
}
\references{
Eddy, W.F. (1982). Convex hull peeling. In: Caussinus, H., Ettinger, P. and Tomassone, R. (eds), \emph{COMPSTAT 1982. Proceedings in computational statistics}, Physica-Verlag (Vienna), 42--47.

Cascos, I. (2009). Data depth: multivariate statistics and geometry. In: Kendall, W.S. and Molchanov, I. (eds) \emph{New Perspectives in Stochastic Geometry}, Clarendon/Oxford University Press (Oxford).
}
\seealso{
\code{\link{depth.halfspace}} for calculation of the Tukey depth.

\code{\link{depth.L2}} for calculation of L2-depth.

\code{\link{depth.Mahalanobis}} for calculation of Mahalanobis depth.

\code{\link{depth.projection}} for calculation of projection depth.

\code{\link{depth.simplicial}} for calculation of simplicial depth.

\code{\link{depth.simplicialVolume}} for calculation of simplicial volume depth.

\code{\link{depth.spatial}} for calculation of spatial depth.

\code{\link{depth.potential}} for calculation of data potential.

\code{\link{depth.zonoid}} for calculation of zonoid depth.

}
\examples{
# Mixture of 3-variate normal distributions
data <- mvrnorm(25, rep(0, 3), diag(3))
x <- rbind(mvrnorm(10, rep(1, 3), diag(3)), data)
depths <- depth.qhpeeling(x, data)
cat("Depths:", depths, "\n")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
