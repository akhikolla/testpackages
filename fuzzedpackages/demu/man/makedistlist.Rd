\name{makedistlist}
\alias{makedistlist}
\title{Make list of distance matrices for calculating GP correlation matrices.}
\usage{
makedistlist(design)
}
\arguments{
\item{design}{An \eqn{n \times p} matrix of input settings}

}
\value{
A list of \eqn{p} matrices, each of dimension \eqn{n \times n} that contain the outer subtractions of each variable in the design matrix.
}
\description{
\code{makedistlist()} is a helper function used to setup the difference matrices that are used by the DPP models. 
}
\examples{
library(demu)

design=matrix(runif(10,1,5),ncol=2,nrow=5)
r=getranges(design)
design=scaledesign(design,r)
l.v=makedistlist(design)
}
\seealso{
  \code{\link{getranges}}
  \code{\link{scaledesign}}
}
