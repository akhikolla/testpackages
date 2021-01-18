\name{ftf}
\alias{ftf}
\title{Fast Trend Filtering}
\description{
  This function solves the structural filtering problem via the primal dual active set algorithm. It fits a k-th order piecewise polynomial by minimizing the least squares error with constraints on the number of breaks in their (k + 1)-st discrete derivative, for a chosen integer k >= 0.
}

\usage{
 ftf(y, k = 1, s = 20, K.max = 5)
}
\arguments{
  \item{y}{Numeric vector of inputs.}
  \item{k}{An integer specifying the desired order of the piecewise polyomial produced by the solution of the trend filtering problem. Must be non-negative, and the default to 1 (linear trend filtering).}
  \item{s}{Number of knots in the piecewise polynomial(breaks in the (k+1)-st derivative), default is 20.}
  \item{K.max}{The maximum number of steps for the algorithm to take before termination. Default is 5.}
}
\details{
  The L0 trend filtering fits an adaptive piecewise polynomial to linearly ordered observations with contraints on the number of knots, for a chosen integer k >= 0. The knots or the breaks in their (k + 1)-st discrete derivative are chosen adaptively based on the observations.

}

\value{
  \item{y}{The observed response vector. Useful for plotting and other methods.}
  \item{beta}{Filtered value.}
  \item{v}{Primal coefficient. The indexes of the nonzero values correspond to the locations of the breaks. }
}



\author{Canhong Wen, Xueqin Wang, Yanhe Shen, Aijun Zhang}
\references{
 Wen,C., Wang, X., Shen, Y., and Zhang, A. (2017). "L0 trend filtering", technical report.
}


\seealso{
\code{\link{plotl0}}.
}


\examples{
set.seed(1)

sigma <- 0.5
y0 <- c((10:30)/3, (40:10)/4, 2:8)
y <- y0 + sigma*rnorm(length(y0))
re <- ftf(y, k = 1, s = 5)


}



