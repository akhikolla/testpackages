\name{fsfused.ada}
\alias{fsfused.ada}
\title{Adaptive Fast Sparse Fused Regression}
\description{
  This is a function that solves the L0 fused problem via the primal dual active set algorithm in sparse condition. It fits a piecewise constant regression model by minimizing the number of breaks in derivative with constraints on the least squares error.}

\usage{
 fsfused.ada(y, tau=1, s.max=20, T, eps=0.1)
}
\arguments{
  \item{y}{Numeric vector of inputs.}
  \item{tau}{Step length for searching the best model, i.e., in the t-th iteration, a model with tau*t knots will be fitted.}
  \item{s.max}{The maximum nubmer of knots in the piecewise constant(breaks in the (k+1)-st derivative), default is 20}
  \item{T}{Number of non-zero values in fitted coefficient.}
  \item{eps}{Early stop criterion. The algorithm stops when mean squared error is less than eps}
}
\value{
  \item{y}{The observed response vector. Useful for plotting and other methods.}
  \item{beta}{Fitted value}
  \item{v}{Primal coefficient. The indexes of the nonzero values correspond to the locations of the breaks. }
  \item{beta.all}{Solution path of fitted value, beta, corresponding to different degrees of freedom.}
  \item{df}{A vector giving an unbiased estimate of the degrees of freedom of the fit, i.e., the number of nonzero values in \code{v}.}
}

\author{Canhong Wen, Xueqin Wang, Yanhe Shen, Aijun Zhang}
\references{
 Wen,C., Wang, X., Shen, Y., and Zhang, A. (2017). "L0 trend filtering", technical report.
}


\seealso{
\code{\link{fsfused}}.
}


\examples{
n <- 1000
sigma <- 0.5
y0 <- rep(0,n)
y0[100:150] <- 2.5
y0[400:600] <- -2.4
y0[800:810] <- 4
y <- y0 + sigma*rnorm(n)

re = fsfused.ada(y, tau=1, s.max=10, T = 260, eps=1.2*sigma^2)


}



