\name{gen.phiCI}
\alias{gen.phiCI}
\title{Average kinship confidence interval}
\description{Returns the confidence interval of the average kinship among pairs of specified individuals}
\usage{gen.phiCI(phiMatrix, prob=c(0.025,0.05,0.95,0.975), b=5000, print.it=F)}
\arguments{
\item{phiMatrix}{A square matrix of kinship coefficients obtained for example with gen.phi, 
	or an array of square matrices (GLmultiMatrix object) of kinship coefficients obtained with gen.phi. Required.}
\item{prob}{Confidence limits probabilities. Default is probs=c(0.025, 0.05, 0.95, 0.975).}
\item{b}{Number of simulations used to calculate the confidence interval. Default is b="5000".}
\item{print.it}{If print.it="F" (default), progression in the number of simulations performed is not displayed.}
}
\value{A vector or a GLmultiVector object depending on the type of 'phiMatrix' entered.
	  If 'phiMatrix' is a GLmultiMatrix, an array of vectors (GLmultiVector), one for each depth, is returned.
	  Array of size ('depthMax'-'depthMin'), vectors of size 'length(probs)'.
	  If 'phiMatrix' is a matrix, a vector of size 'length(probs)' is returned.}
\seealso{

\code{\link{gen.genealogy}}
\code{\link{gen.phi}}
\code{\link{gen.phiOver}}
\code{\link{gen.phiMean}}
}
\examples{
\donttest{
data(genea140)
gen140<-gen.genealogy(genea140)
ge07<-gen.branching(gen140, pro=gen.pro(gen140)[c(1:7)])
phiMat<-gen.phi(ge07)
gen.phiCI(phiMat)
}
}
\keyword{manip}
