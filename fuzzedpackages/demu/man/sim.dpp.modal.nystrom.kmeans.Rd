\name{sim.dpp.modal.nystrom.kmeans}
\alias{sim.dpp.modal.nystrom.kmeans}
\title{Subsample an observational dataset using the conditional DPP design emulator with a kmeans-based Nystrom approximation.}
\usage{
sim.dpp.modal.nystrom.kmeans(Xin,rho=rep(0.01,ncol(Xin)),
  n,m=max(ceiling(nrow(Xin)*0.1),n),method="KmeansNystrom",
  initializer="kmeans++",...)
}
\arguments{
\item{Xin}{An \eqn{n\times p} dataset of observations from which we want to draw subsamples.}

\item{n}{Size of the designed subsample to draw from \code{Xall}.}

\item{rho}{The \eqn{p\times 1} parameter vector for the Gaussian correlation model.}

\item{m}{Number of landmark points to use in constructing the kmeans-based Nystrom approximation.}

\item{method}{Type of approximation to use.  Currently only supports ``KmeansNystrom''.}

\item{initializer}{Initialization to use in the Kmeans algorithm, default is ``kmeans++''.}

\item{...}{Additional options to pass to \code{MiniBatchKmeans()} for selecting the landmark points.}
}
\value{
A list containing a matrix which is the union of the observation matrix \code{Xin} and selected landmark sites, the indices into this matrix of the selected design sites as well as matrix of the design sites.
}
\description{
\code{sim.dpp.modal.nystrom.kmeans()} uses the kmeans-based Nystrom approximation of Zhang and Kwok (2010) to select \code{n} design sites from the observational dataset \code{Xin} using the DPP-based design emulator of Pratola et al. (2018).  The design constructed assumes a Gaussian process
regression model with stationary correlation function \eqn{r(x,x^\prime)}, where the
entries of \code{R} are formed by evaluating \eqn{r(x,x^\prime)} over a set of landmarks chosen by the kmeans algorithm, and the resulting eigenvectors are projected into the higher dimensional space using the Nystrom approximation.  Additional options for the \code{MiniBatchKmeans()} algorithm from package \code{ClusterR} can be passed to alter the construction of the landmark set.
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

# Fake dataset in 5 dimensions
X=matrix(runif(500*5),ncol=5)
rho=rep(0.01,5)
n=50
samp=sim.dpp.modal.nystrom.kmeans(X,rho,n)
samp$design

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
  \code{\link{sim.dpp.modal.nystrom}}
}
