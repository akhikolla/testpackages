\name{tensorCentering}
\alias{tensorCentering}
\title{Center an Array of Observations}
\description{
Centers an array of array-valued observations by substracting a location array (the mean array by default) from each observation.
}
\usage{
tensorCentering(x, location = NULL)
}
\arguments{
  \item{x}{Array of order at least two with the last dimension corresponding to the sampling units.}
  \item{location}{The location to be used in the centering. Either \code{NULL}, defaulting to the mean array, or a user-specified \eqn{p_1 \times p_2 \times \ldots \times p_r}{p_1 x p_2 x ... x p_r}-dimensional array.}
}
\details{
Centers a \eqn{p_1 \times p_2 \times \ldots \times p_r \times n}{p_1 x p_2 x ... x p_r x n}-dimensional array by substracting the \eqn{p_1 \times p_2 \times \ldots \times p_r}{p_1 x p_2 x ... x p_r}-dimensional \code{location} from each of the observed arrays.
}
\value{
Array of centered observations with the same dimensions as the input array. The used location is returned as attribute \code{"location"}.
}
\author{Joni Virta}
\examples{
## Generate sample data.
n <- 1000
x <- t(cbind(rnorm(n, mean = 0),
             rnorm(n, mean = 1),
             rnorm(n, mean = 2),
             rnorm(n, mean = 3),
             rnorm(n, mean = 4),
             rnorm(n, mean = 5)))
             
dim(x) <- c(3, 2, n)

## Centered data
xcen <- tensorCentering(x)

## Check the means of individual cells
apply(xcen, 1:2, mean)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ array }
