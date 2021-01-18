\name{ffused}
\alias{ffused}
\title{Fast Fused Regression}
\description{
    This is a fast calculation function that solves the L0 fused problem via the primal dual active set algorithm. It fits a piecewise constant regression model by minimizing the least squares error with constraints on the number of breaks in their 1-st discrete derivative.
}

\usage{
 ffused(y, s, K.max=5)
}
\arguments{
  \item{y}{Response sequence to be fitted.}
  \item{s}{Number of knots in the piecewise constant(breaks in the derivative), default is 10.}
  \item{K.max}{The maximum number of steps for the algorithm to take before termination. Default is 5.}

}
\value{
  \item{y}{The observed response vector. Useful for plotting and other methods.}
  \item{beta}{Fitted value.}
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

set.seed(111)
n <- 1000
sigma <- 0.5
y0 <- rep(0,n)
y0[100:150] <- 2.5
y0[400:600] <- -2.4
y0[800:810] <- 4
y <- y0 + sigma*rnorm(n)

re = ffused(y, s = 8, K.max = 5)


}

