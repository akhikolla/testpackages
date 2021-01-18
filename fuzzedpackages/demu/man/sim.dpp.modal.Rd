\name{sim.dpp.modal}
\alias{sim.dpp.modal}
\title{Draw samples from the conditional DPP design emulator.}
\usage{
sim.dpp.modal(R,n=0,eigs=NULL)
}
\arguments{
\item{R}{A correlation matrix evaluated over a grid of candidate design sites.}

\item{n}{Size of the design to sample.}

\item{eigs}{One can alternatively pass the pre-computed eigendecomposition of the correlation matrix instead of \code{R}.}

}
\value{
A vector of indices to the sampled design sites.
}
\description{
\code{sim.dpp.modal()} uses the DPP-based design emulator of Pratola et al. (2018)
to draw a sample of the \code{n}-run optimal design for a Gaussian process
regression model with stationary correlation function \eqn{r(x,x^\prime)}, where the
entries of \code{R} are formed by evaluating \eqn{r(x,x^\prime)} over a grid of candidate 
locations.
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

# candidate grid
ngrid=20
x=seq(0,1,length=ngrid)
X=as.matrix(expand.grid(x,x))
l.d=makedistlist(X)

# draw design from DPP mode
n=21
rho=0.01
R=rhomat(l.d,rep(rho,2))$R
pts=sim.dpp.modal(R,n)

# Could plot the result:
# plot(X,xlim=c(0,1),ylim=c(0,1))
# points(X[pts,],pch=20)
}
\seealso{
  \code{\link{demu-package}}
  \code{\link{sim.dpp.modal.fast}}
  \code{\link{sim.dpp.modal.seq}}
}
