\name{fsf.ada}
\alias{fsf.ada}
\title{Adaptive Fast Structural Filtering}
\description{
  This is a function that solves the structural filtering problem with L0 penalty via the primal dual active set algorithm. It fits a non-parametric regression model by minimizing the number of nonzero values in D*beta with constraints on the least squares error.
}

\usage{
 fsf.ada(y, D, tau=1, s.max=20,  eps=0.1, ddinv=NULL)
}
\arguments{
  \item{y}{Response sequence to be filtered.}
  \item{D}{Penalty matrix on coeffient beta.}
  \item{tau}{Step length for searching the best model, i.e., in the t-th iteration, a model with tau*t knots will be fitted.}
  \item{s.max}{The maximum nubmer of knots in the penalized coefficient(breaks in \code{D*beta}, default is 20}
  \item{eps}{Early stop criterion. The algorithm stops when mean squared error is less than eps.}
  \item{ddinv}{The inverse matrix of \code{D*t(D)}, could be NULL input.}
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
ddt <- D0\%*\%t(D0)
ddinv<- Solve.banded(ddt, 1,1, B = diag(1,dim(D0)[1]))

re <- fsf.ada(y, D0, tau = 1, s.max = 10, eps = 0.1, ddinv = ddinv)

}


