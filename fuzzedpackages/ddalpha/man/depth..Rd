\name{depth.}
\alias{depth.}
\title{
Calculate Depth
}
\description{
Calculates the depth of points w.r.t. a multivariate data set.

The detailed descriptions are found in the corresponding topics.
}
\usage{
depth.(x, data, notion, ...)

## beta-skeleton depth
# depth.betaSkeleton(x, data, beta = 2, distance = "Lp", Lp.p = 2, 
#                   mah.estimate = "moment", mah.parMcd = 0.75)

## Tukey depth
# depth.halfspace(x, data, exact, method, num.directions = 1000, seed = 0)

## L2-depth
# depth.L2(x, data, mah.estimate = "moment", mah.parMcd = 0.75)

## Mahalanobis depth
# depth.Mahalanobis(x, data, mah.estimate = "moment", mah.parMcd = 0.75)

## projection depth
# depth.projection(x, data, method = "random", num.directions = 1000)

## simplicial depth
# depth.simplicial(x, data, exact = F, k = 0.05, seed = 0)

## simplicial volume depth
# depth.simplicialVolume(x, data, exact = F, k = 0.05, seed = 0)

## spatial depth
# depth.spatial(x, data)

## zonoid depth
# depth.zonoid(x, data)

## potential
# depth.potential (x, data, pretransform = "1Mom", 
#            kernel = "GKernel", kernel.bandwidth = NULL, mah.parMcd = 0.75)

## convex hull peeling depth
# depth.qhpeeling(x, data)

}
\arguments{
  \item{x}{
Matrix of objects (numerical vector as one object) whose depth is to be calculated; each row contains a \eqn{d}-variate point. Should have the same dimension as \code{data}.
}
  \item{data}{
Matrix of data where each row contains a \eqn{d}-variate point, w.r.t. which the depth is to be calculated.
}
  \item{notion}{
The name of the depth notion (shall also work with a user-defined depth function named \code{"depth.<name>"}).
}
  \item{\dots}{
Additional parameters passed to the depth functions.
}
}

\seealso{

\code{\link{depth.betaSkeleton}}

\code{\link{depth.halfspace}}

\code{\link{depth.L2}}

\code{\link{depth.Mahalanobis}}

\code{\link{depth.projection}}

\code{\link{depth.simplicial}}

\code{\link{depth.simplicialVolume}}

\code{\link{depth.spatial}}

\code{\link{depth.zonoid}}

\code{\link{depth.potential}}

\code{\link{depth.qhpeeling}}

\code{\link{depth.graph}} for building the depth surfaces of the two dimensional data.

}
\value{
Numerical vector of depths, one for each row in \code{x}; or one depth value if \code{x} is a numerical vector.
}
\examples{
# 5-dimensional normal distribution
data <- mvrnorm(1000, rep(0, 5), 
                matrix(c(1, 0, 0, 0, 0, 
                         0, 2, 0, 0, 0, 
                         0, 0, 3, 0, 0, 
                         0, 0, 0, 2, 0, 
                         0, 0, 0, 0, 1),
                nrow = 5))
x <- mvrnorm(10, rep(1, 5), 
             matrix(c(1, 0, 0, 0, 0, 
                      0, 1, 0, 0, 0, 
                      0, 0, 1, 0, 0, 
                      0, 0, 0, 1, 0, 
                      0, 0, 0, 0, 1),
             nrow = 5))
                
depths <- depth.(x, data, notion = "zonoid")
cat("Depths: ", depths, "\n")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
