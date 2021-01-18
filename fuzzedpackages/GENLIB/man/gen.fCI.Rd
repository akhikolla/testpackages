\name{gen.fCI}
\alias{gen.fCI}
\title{Average inbreeding coefficient confidence interval}
\description{Returns the confidence interval of the average inbreeding among specified individuals}
\usage{gen.fCI(vectF, prob=c(0.025,0.05,0.95,0.975), b="5000", print.it="F")}
\arguments{
\item{vectF}{A vector of inbreeding coefficients obtained for example with gen.f. Required.}
\item{prob}{Confidence limits probabilities. Default is probs=c(0.025, 0.05, 0.95, 0.975).}
\item{b}{Number of simulations used to calculate the confidence interval. Default is b="5000".}
\item{print.it}{If print.it="F" (default), progression in the number of simulations performed is not displayed.}
}
\value{A vector or a GLmultiVector object depending on the type of 'vectF' entered.
	  If 'vectF' is a GLmultiVector, an array of vectors (GLmultiVector), one for each depth, is returned.
	  Array of size ('depthMax'-'depthMin'), vectors of size 'length(probs)'.
	  If 'vectF' is a GLmultiVector, a vector of size 'length(probs)' is returned.}
\seealso{
\code{\link{gen.f}}
}
\examples{
data(genea140) 
gen140<-gen.genealogy(genea140) 
gen10 <-gen.branching(gen140, pro=gen.pro(gen140)[c(1:10)])
fval<-gen.f(gen10)
gen.fCI(fval)
}
\keyword{manip}
