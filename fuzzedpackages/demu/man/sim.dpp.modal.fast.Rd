\name{sim.dpp.modal.fast}
\alias{sim.dpp.modal.fast}
\title{Draw samples from the conditional DPP design emulator.}
\usage{
sim.dpp.modal.fast(R,n)
}
\arguments{
\item{R}{A sparse correlation matrix evaluated over a grid of candidate design sites.  The sparse matrix should be of type \code{dgCMatrix} (see package \link{spam}).}

\item{n}{Size of the design to sample.}


}
\value{
A vector of indices to the sampled design sites.
}
\description{
\code{sim.dpp.modal.fast()} is similar to \code{sim.dpp.modal} but is a C++ codepath that makes use of
SPAM's sparse matrices to enable faster computation.  It implements the DPP-based design 
emulator of Pratola et al. (2018) to draw a sample of the \code{n}-run optimal design for a Gaussian 
process regression model with compact correlation function \eqn{r(x,x^\prime)}, where the entries 
of \code{R} are formed by evaluating \eqn{r(x,x^\prime)} over a grid of candidate locations.
}
\details{
For more details on the method, see Pratola et al. (2018).  Detailed
examples demonstrating the method are available at
\url{http://www.matthewpratola.com/software}.
}
\references{
  Pratola, Matthew T., Lin, C. Devon, and Craigmile, Peter. (2018)
  Optimal Design Emulators: A Point Process Approach. 
  \emph{arXiv:1804.02089}.
}
\examples{
library(demu)
library(fields)
library(spam)
library(Matrix)
library(Rcpp)

# candidate grid
ngrid=20
x=seq(0,1,length=ngrid)
X=as.matrix(expand.grid(x,x))

# draw design from DPP mode
n=21
theta=0.39
R.spam=wendland.cov(X,X,theta=theta,k=3) 
R=as.dgCMatrix.spam(R.spam)
rm(R.spam)
pts=sim.dpp.modal.fast(R,n)

# Could plot the result:
# plot(X,xlim=c(0,1),ylim=c(0,1))
# points(X[pts,],pch=20)
}
\seealso{
  \code{\link{demu-package}}
  \code{\link{sim.dpp.modal}}
  \code{\link{sim.dpp.modal.seq}}
}
