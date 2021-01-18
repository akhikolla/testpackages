\name{sim.dpp.modal.np}
\alias{sim.dpp.modal.np}
\title{Draw samples from the conditional DPP design emulator using a kmeans-based Nystrom approximation.}
\usage{
sim.dpp.modal.np(n,p,N,rho,m=max(ceiling(N*0.1),n),...)
}
\arguments{
\item{n}{Size of the desired design.}

\item{p}{Dimension of the desired design.}

\item{N}{Number of kernel approximation points drawn uniformly from the \code{p}-dimensional design space.}

\item{rho}{The \eqn{p\times 1} parameter vector for the Gaussian correlation model.}

\item{m}{Number of landmark points to use in constructing the kmeans-based Nystrom approximation.}

\item{...}{Additional options to pass to \code{sim.dpp.modal.nystrom.kmeans()} for drawing the design.}

}
\value{
A list containing a matrix which is the union of the \eqn{N\times p} uniformly sampled kernel approximation points and the \code{m} selected landmark sites, and the indices into this matrix of the selected design sites.
}
\description{
\code{sim.dpp.modal.np()} uses \code{sim.dpp.modal.nystrom.kmeans()} to draw a design of \code{n} points in \code{p} dimensions using the kmeans-based Nystrom approximation of Zhang and Kwok (2010) and the DPP-based design emulator of Pratola et al. (2018).  The design constructed assumes a Gaussian process
regression model with stationary correlation function \eqn{r(x,x^\prime)}, where the
entries of \code{R} are formed by evaluating \eqn{r(x,x^\prime)} over a set of landmarks chosen by the kmeans algorithm, and the resulting eigenvectors are projected into the higher dimensional space using the Nystrom approximation.  Additional options for \code{sim.dpp.modal.nystrom.kmeans()} can be passed to alter the construction of the landmark set.
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

  Zhang, Kai and Kwok, James T. (2010)
  Clustered Nystrom method for large scale manifold learning and dimension reduction.
  \emph{IEEE Transactions on Neural Networks}, \bold{21.10}, 1576--1587.
  \doi{10.1109/TNN.2010.2064786}
}
\examples{
library(demu)

n=50
p=5
N=500
rho=rep(0.01,5)
samp=sim.dpp.modal.np(n,p,N,rho)

# Could plot the result:
# pchvec=rep(1,nrow(samp$X))
# pchvec[samp$pts]=20
# cexvec=rep(0.1,nrow(samp$X))
# cexvec[samp$pts]=1
# colvec=rep("black",nrow(samp$X))
# colvec[samp$pts]="red"
# pairs(samp$X,pch=pchvec,cex=cexvec,col=colvec,xlim=c(0,1),ylim=c(0,1))
}
\seealso{
  \code{\link{demu-package}}
  \code{\link{sim.dpp.modal}}
  \code{\link{sim.dpp.modal.nystrom.kmeans}}
}
