\name{sim.dpp.modal.fast.seq}
\alias{sim.dpp.modal.fast.seq}
\title{Draw sequential samples from the conditional DPP given previously sampled points already in the design.}
\usage{
sim.dpp.modal.fast.seq(curpts, R,n)
}
\arguments{
\item{curpts}{A vector of indices to the candidate points that already appear in the design.}

\item{R}{A sparse correlation matrix evaluated over a grid of candidate design sites.  The sparse matrix should be of type \code{dgCMatrix} (see package \link{spam}).}

\item{n}{Size of the design to sample.}


}
\value{
A vector of indices to the sampled design sites.
}
\description{
\code{sim.dpp.modal.fast.seq()} is similar to \code{sim.dpp.modal.fast} but sequentially selects \code{n} additional points to 
add to the design given that the points in \code{curpts} are alread in the design from previous sequential 
iterations.  
It uses the C++ codepath that makes use of
SPAM's sparse matrices to enable faster computation.  It implements the DPP-based design 
emulator of Pratola et al. (2018) to draw a sequential sample of the \code{n}-run additional optimal design points for a Gaussian 
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

n1=3
n2=3
n3=3
rho=0.2
ngrid=10

x=seq(0,1,length=ngrid)
X=as.matrix(expand.grid(x,x))
l.d=makedistlist(X)

# Initial design
R.spam=wendland.cov(X,X,theta=rho,k=3)
R=as.dgCMatrix.spam(R.spam)
pts.1=sim.dpp.modal.fast(R,n1)
pts.1.proj=remove.projections(pts.1,X)

# Next sequential step, removing projections
pts.2=sim.dpp.modal.fast.seq(pts.1.proj$allpts,R,n2)
design=c(pts.1,pts.2$pts.new)
pts.2.proj=remove.projections(design,X)

# Next sequential step, removing projections
pts.3=sim.dpp.modal.fast.seq(pts.2.proj$allpts,R,n3)
design=c(design,pts.3$pts.new)


# Or, starting with the initial design, don't remove projections
pts.2=sim.dpp.modal.fast.seq(pts.1,R,n2)
designB=c(pts.1,pts.2$pts.new)

pts.3=sim.dpp.modal.fast.seq(designB,R,n3)
designB=c(designB,pts.3$pts.new)


# Plot the result:
#par(mfrow=c(1,3))
#plot(X,xlim=c(0,1),ylim=c(0,1),main="Initial Design")
#points(X[pts.1,],pch=20,cex=2)
#
#plot(X,xlim=c(0,1),ylim=c(0,1),main="+3x2 remove projections")
#points(X[design,],pch=20,cex=2)
#
#plot(X,xlim=c(0,1),ylim=c(0,1),main="+3x2 not removing projections")
#points(X[designB,],pch=20,cex=2)
}
\seealso{
  \code{\link{demu-package}}
  \code{\link{sim.dpp.modal.fast}}
  \code{\link{sim.dpp.modal}}
}
