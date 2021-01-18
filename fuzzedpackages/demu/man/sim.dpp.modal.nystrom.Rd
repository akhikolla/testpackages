\name{sim.dpp.modal.nystrom}
\alias{sim.dpp.modal.nystrom}
\title{Draw samples from the conditional DPP design emulator using grid-based Nystrom approximation.}
\usage{
sim.dpp.modal.nystrom(Xin,rho,n=0,ngrid=NULL,method="Nystrom")
}
\arguments{
\item{Xin}{A initial \eqn{n\times p} matrix of points.}

\item{rho}{The \eqn{p\times 1} parameter vector for the Gaussian correlation model.}

\item{n}{Size of the design to sample from the candidate grid.}

\item{ngrid}{Size of the candidate grid will be \code{ngrid^p}.}

\item{method}{Type of approximation to use.  Currently only supports ``Nystrom''.}

}
\value{
A list containing the candidate points constructed and the points selected as the design sites from this candidate set as well as their indices.
}
\description{
\code{sim.dpp.modal.nystrom()} uses the DPP-based design emulator of Pratola et al. (2018)
to draw a sample of the \code{n}-run optimal design for a Gaussian process
regression model with stationary correlation function \eqn{r(x,x^\prime)}, where the
entries of \code{R} are formed by evaluating \eqn{r(x,x^\prime)} over a grid of candidate 
locations.  This function uses a grid-based Nystrom approximation based on the passed matrix \code{X} to avoid constructing a large correlation matrix if dimension \code{ngrid^p} and its subsequent eigendecomposition.
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

# Starting design
X=matrix(runif(10*2),ncol=2)
rho=rep(0.01,2)
n=10
ngrid=11
samp=sim.dpp.modal.nystrom(X,rho,n,ngrid)
samp$design

# Could plot the result:
# plot(samp$X,xlim=c(0,1),ylim=c(0,1))
# points(samp$X[samp$pts,],pch=20)
}
\seealso{
  \code{\link{demu-package}}
  \code{\link{sim.dpp.modal}}
  \code{\link{sim.dpp.modal.nystrom.kmeans}}
}
