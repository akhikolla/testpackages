\name{plot.FastRCS}
\alias{plot.FastRCS}
\title{Robust Diagnostic Plots For FastRCS}
\description{
  Shows the robust Score distances versus robust Orthogonal distances
  and their respective cutoffs, for the an object of class FastRCS. 

}
\usage{
\method{plot}{FastRCS}(x,col="black",pch=16,...)
}
\arguments{
  \item{x}{For the \code{plot()} method, a \code{FastRCS} object, typically
    result of \code{\link{FastRCS}}.}	
	\item{col}{ A specification for the default plotting color. Vector of values are recycled.}
	\item{pch}{ Either an integer specifying a symbol or a single character
          to be used as the default in plotting points.   Note that only
          integers and single-character strings can be set as a
          graphics parameter.
          Vector of values are recycled.}
   \item{...}{Further arguments passed to the plot function.}
}
\details{
  This function produces the robust standardized, residuals
	 as well as an indicative cut-off (under normal model). 
	This tool is a diagnostic plot for robust regression and 
	can be used used to reveal the outliers. 
}
\seealso{
  \code{\link{FastRCS}}
}
\examples{
set.seed(123)
n<-100
p<-3
x0<-matrix(rnorm(n*p),nc=p)
y0<-rnorm(n)
z<-c(rep(0,30),rep(1,70))
x0[1:30,]<-matrix(rnorm(30*p,5,1/100),nc=p)
y0[1:30]<-rnorm(30,10,1/100)
ns<-FRCSnumStarts(p=p,eps=0.4);
results<-FastRCS(x=x0,y=y0,alpha=0.5,nSamp=ns)
plot(results)

}
\keyword{hplot}
\keyword{robust}
\keyword{regression}
