\name{plot.FastPCS}
\alias{plot.FastPCS}
\title{Distance plot for FastPCS}
\description{
  Plots the robust distance values from a FastPCS model fit, and their parametric cut-off. 

}
\usage{
\method{plot}{FastPCS}(x,col="black",pch=16,...)
}
\arguments{
  \item{x}{For the \code{plot()} method, a \code{FastPCS} object, typically
    resulting as output from \code{\link{FastPCS}}.}	
	\item{col}{ A specification for the default plotting color. Vectors of values are recycled.}
	\item{pch}{ Either an integer specifying a symbol, or a single character
          to be used as the default in plotting points. Note that only
          integers and single-character strings can be set as
          graphics parameters.
          Vectors of values are recycled.}
   \item{...}{Further arguments passed to the plot function.}
}

\seealso{
  \code{\link{FastPCS}}
}
\examples{
## generate data
set.seed(123)
n<-100
p<-3
x0<-matrix(rnorm(n*p),nc=p)
x0[1:30,]<-matrix(rnorm(30*p,4.5,1/100),nc=p)
z<-c(rep(0,30),rep(1,70))
nstart<-FPCSnumStarts(p=p,eps=0.4)
results<-FastPCS(x=x0,nSamp=nstart)
colvec<-rep("orange",length(z))
colvec[z==1]<-"blue"
plot.FastPCS(results,col=colvec,pch=16)
}
\keyword{hplot}
\keyword{robust}
\keyword{multivariate}
\author{Kaveh Vakili, Eric Schmitt}
