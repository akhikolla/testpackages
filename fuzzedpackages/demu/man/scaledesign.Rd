\name{scaledesign}
\alias{scaledesign}
\title{Rescale a design matrix to the [0,1] hypercube.}
\usage{
scaledesign(design,r)
}
\arguments{
\item{design}{An \eqn{n \times p} matrix of input settings}
\item{r}{An \eqn{p \times 2} matrix of variable ranges extracted from \code{getranges()}}

}
\value{
A \eqn{n \times p} design matrix rescaled to the \eqn{[0,1]} hypercube.
}
\description{
\code{scaledesign()} is a helper function to rescale a design to the \eqn{[0,1]} hypercube using variable ranges previously extracted by a call to \code{getranges()}.
}
\examples{
library(demu)

design=matrix(runif(10,1,5),ncol=2,nrow=5)
r=getranges(design)
scaledesign(design,r)

}
\seealso{
  \code{\link{unscalemat}}
}
