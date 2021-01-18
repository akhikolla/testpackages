\name{fsf}
\alias{fsf}
\title{Fast Structural Filtering}
\description{
  This function solves the generalized structural filtering problem via the primal dual active set algorithm. It fits a non-parametric regression model by minimizing the least squares error with penalty matrix D on coefficient beta.
}

\usage{
 fsf(y, D, s = 20, K.max = 5, ddinv=NULL)
}
\arguments{
  \item{y}{Response sequence to be filtered.}
  \item{D}{Penalty matrix on coeffient beta.}
  \item{s}{Number of knots in the penalized coefficient(breaks in the \code{D*beta}), default is 20.}
  \item{K.max}{The maximum number of steps for the algorithm to take before termination. Default is 5.}
  \item{ddinv}{The inverse matrix of \code{D*t(D)}, could be NULL input.}
}
\value{
   \item{y}{The observed response vector. Useful for plotting and other methods.}
  \item{beta}{Fitted value.}
  \item{v}{Primal coefficient. The indexes of the nonzero values correspond to the locations of the breaks in \code{D*beta}. }
}


\author{Canhong Wen, Xueqin Wang, Yanhe Shen, Aijun Zhang}
\references{
 Wen,C., Wang, X., Shen, Y., and Zhang, A. (2017). "L0 trend filtering", technical report.
}


\seealso{
\code{\link{plotl0}}.
}



\examples{
require(limSolve)
n <- 1000
sigma <- 0.5
y0 <- rep(0,n)
y0[100:150] <- 2
y0[400:600] <- -1
y0[800:810] <- 4
y <- y0 + sigma*rnorm(n)
y[800:810] <- y0[800:810] + sigma*rnorm(11)
D0 <- matrix(0, n-1,n)
diag(D0) <- -1
for(i in 1:(n-1)) D0[i,i+1] <- 1
m <- dim(D0)[1]

re = fsf(y, D0)

}




