\name{unscalemat}
\alias{unscalemat}
\title{Unscale a matrix back to its original ranges.}
\usage{
unscalemat(mat,r)
}
\arguments{
\item{mat}{An \eqn{n \times p} matrix of numbers scaled to the \eqn{[0,1]} hypercube}
\item{r}{An \eqn{p \times 2} matrix of the original ranges of the variables}

}
\value{
A \eqn{n \times p} matrix with variables rescaled back to their original ranges, as specified by \code{ranges}.
}
\description{
\code{unscalemat()} is a helper function to rescale a matrix back to its original ranges.  Typically this is used to rescale the posterior samples of the parameters back to their original scale.
}
\examples{
library(demu)

design=matrix(runif(10,1,5),ncol=2,nrow=5)
r=getranges(design)
design=scaledesign(design,r)
unscalemat(design,r)

}
\seealso{
  \code{\link{getranges}}
  \code{\link{scaledesign}}
}
