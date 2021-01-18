\name{ftf.ada}
\alias{ftf.ada}
\title{Adaptive Fast Trend Filtering}
\description{
  This is a function that adaptively solves the trend filtering problem with L0 penalty via the primal dual active set algorithm. It fits a k-th order piecewise polynomial by minimizing the number of breaks in the (k + 1)-st discrete derivative with the constraints on the least squares error.
}

\usage{
 ftf.ada(y, k = 1, tau = 1, s.max=20, eps=0.1)
}
\arguments{
  \item{y}{Numeric vector of inputs.}
  \item{k}{An integer specifying the desired order of the piecewise polyomial produced by the solution of the trend filtering problem. Must be non-negative, and the default to 1 (linear trend filtering).}
  \item{tau}{Step length for searching the best model, i.e., in the t-th iteration, a model with tau*t knots will be fitted.}
  \item{s.max}{The maximum nubmer of knots in the piecewise polynomial(breaks in the (k+1)-st derivative), default is 20}
  \item{eps}{Early stop criterion. The algorithm stops when mean squared error is less than eps}
}
\details{
  The L0 trend filtering fits an adaptive piecewise polynomial to linearly ordered observations with contraints on the number of knots, for a chosen integer k >= 0. The knots or the breaks in their (k + 1)-st discrete derivative are chosen adaptively based on the observations.

}

\value{
  \item{y}{The observed response vector. Useful for plotting and other methods.}
  \item{beta}{Filtered value}
  \item{v}{Primal coefficient. The indexes of the nonzero values correspond to the locations of the breaks. }
  \item{beta.all}{Solution path of filtered value, beta, corresponding to different degrees of freedom.}
  \item{df}{A vector giving an unbiased estimate of the degrees of freedom of the fit, i.e., the number of nonzero values in \code{v}.}
}



\author{Canhong Wen, Xueqin Wang, Yanhe Shen, Aijun Zhang}
\references{
 Wen,C., Wang, X., Shen, Y., and Zhang, A. (2017). "L0 trend filtering", technical report.
}


\seealso{
\code{\link{ftf}}.
}


\examples{
set.seed(1)

sigma <- 0.5
y0 <- c((10:30)/3, (40:10)/4, 2:8)
y <- y0 + sigma*rnorm(length(y0))
re <- ftf.ada(y, k = 1, s.max = 5)


}



