\name{depth.betaSkeleton}
\alias{depth.betaSkeleton}
\title{
Calculate Beta-Skeleton Depth
}
\description{
Calculates the beta-skeleton depth of points w.r.t. a multivariate data set.
}
\usage{
depth.betaSkeleton(x, data, beta = 2, distance = "Lp", Lp.p = 2, 
                   mah.estimate = "moment", mah.parMcd = 0.75)
}
\arguments{
  \item{x}{Matrix of objects (numerical vector as one object) whose depth is to be calculated; each row contains a \eqn{d}-variate point. Should have the same dimension as \code{data}.
}
  \item{data}{Matrix of data where each row contains a \eqn{d}-variate point, w.r.t. which the depth is to be calculated.
}
  \item{beta}{The paremeter defining the positionning of the balls' centers, see Yang and Modarres (2017) for details. By default (together with other arguments) equals \code{2}, which corresponds to the lens depth, see Liu and Modarres (2011).
}
  \item{distance}{A character string defining the distance to be used for determining inclusion of a point into the lens (influence region), see Yang and Modarres (2017) for details. Possibilities are \code{"Lp"} for the Lp-metric (default) or \code{"Mahalanobis"} for the Mahalanobis distance adjustment.
  }
  \item{Lp.p}{A non-negative number defining the distance's power equal \code{2} by default (Euclidean distance); is used only when \code{distance = "Lp"}.
  }
  \item{mah.estimate}{A character string specifying which estimates to use when calculating sample covariance matrix; can be \code{"none"}, \code{"moment"} or \code{"MCD"}, determining whether traditional moment or Minimum Covariance Determinant (MCD) (see \code{\link{covMcd}}) estimates for mean and covariance are used. By default \code{"moment"} is used. Is used only when \code{distance = "Mahalanobis"}.
}
  \item{mah.parMcd}{The value of the argument \code{alpha} for the function \code{\link{covMcd}}; is used when \code{distance = "Mahalanobis"} and \code{mah.estimate = "MCD"}.
}
}
\details{
Calculates the beta-skeleton depth, see Yang and Modarres (2017). Its particular case, lens depth, see Liu and Modarres (2011), is obtained when \code{beta = 2}, \code{distance = "Lp"} and \code{Lp.p = 2} (default settings). For tne example of the lens depth, the depth of an observation \code{x} is calculated as the portion of lens containing \code{x}, with lens being an intersection of two closed balls centered at two sample's points each having radius equal to the distance between these two points.
}
\value{
Numerical vector of depths, one for each row in \code{x}; or one depth value if \code{x} is a numerical vector.
}
\references{
Liu, Z. and Modarres, R. (2011). Lens data depth and median. \emph{Journal of Nonparametric Statistics} \bold{23}(4) 1063--1074.

Yang, M. and Modarres, R. (2017). \eqn{\beta}-skeleton depth functions and medians. \emph{Commmunications in Statistics - Theory and Methods} to appear.
}
\seealso{
\code{\link{depth.halfspace}} for calculation of the Tukey depth.

\code{\link{depth.Mahalanobis}} for calculation of Mahalanobis depth.

\code{\link{depth.projection}} for calculation of projection depth.

\code{\link{depth.simplicial}} for calculation of simplicial depth.

\code{\link{depth.simplicialVolume}} for calculation of simplicial volume depth.

\code{\link{depth.spatial}} for calculation of spatial depth.

\code{\link{depth.zonoid}} for calculation of zonoid depth.

\code{\link{depth.potential}} for calculation of data potential.

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
                
depths <- depth.betaSkeleton(x, data)
cat("Depths:", depths, "\n")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
