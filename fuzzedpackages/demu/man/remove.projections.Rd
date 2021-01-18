\name{remove.projections}
\alias{remove.projections}
\title{Identify candidate points making up all marginal subprojections of an existing design.}
\usage{
remove.projections(curpts,X)
}
\arguments{
\item{curpts}{Indices of points currently in the design}

\item{X}{An \code{n x p} matrix of all candidate points}

}
\value{
A list containing the vector \code{curpts}, the vector \code{projpts} which contains the identified projection points of the current design, and \code{allpts}.
}
\description{
\code{remove.projections()} is a helper function to identify all lower-dimensional marginal projection points of the existing design points indexed by \code{curpts}.  This function can be used to remove a subset of points from the \emph{candidate set} in order to enforce non-collapsingness of when sequentially adding design points.
}
\examples{
library(demu)

n1=3
n2=3
n3=3
rho=rep(1e-10,2)
ngrid=10

x=seq(0,1,length=ngrid)
X=as.matrix(expand.grid(x,x))
l.d=makedistlist(X)

# Initial design
R=rhomat(l.d,rho)$R
pts.1=sim.dpp.modal(R,n1)
pts.1.proj=remove.projections(pts.1,X)

# Plot - design points in black, design+projection points in grey.
#plot(X,xlim=c(0,1),ylim=c(0,1))
#points(X[pts.1.proj$projpts,],pch=20,cex=2,col="grey")
#points(X[pts.1,],pch=20,cex=2)
}
\seealso{
  \code{\link{demu-package}}
  \code{\link{sim.dpp.modal.seq}}
}
