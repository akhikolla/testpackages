\name{plot.FastHCS}
\alias{plot.FastHCS}
\title{Robust diagnostic plots for FastHCS}
\description{
  Creates a diagnostic plot of the robust SD and OD values from a FastHCS model fit, and their parametric cutoffs. 

}
\usage{
\method{plot}{FastHCS}(x,col="black",pch=16,...)
}
\arguments{
  \item{x}{For the \code{plot()} method, a \code{FastHCS} object, typically
    resulting as output from \code{\link{FastHCS}}.}	
	\item{col}{ A specification for the default plotting color. Vectors of values are recycled.}
	\item{pch}{ Either an integer specifying a symbol, or a single character
          to be used as the default in plotting points. Note that only
          integers and single-character strings can be set as
          graphics parameters.
          Vectors of values are recycled.}
   \item{...}{Further arguments passed to the plot function.}
}
\details{
  This function produces the PCA diagnostic plot of Hubert et al. (2005).
 Score distances are the n-vector of distances of each observation 
 to the robust estimate of location on the robust PCA subspace. 
Likewise, orthogonal distances are the n-vector of distances of 
 each observations to the robust PCA subspace. The observations whose
  score distance is larger than cutoff.sd or whose orthogonal distance
 is larger than cutoff.od are considered outliers and receive a
 flag equal to zero.
   The orthogonal distances are displayed along the vertical axis 
and the score distances along the horizontal axis, with the dotted 
 lines indicating their respective cut-offs. 
}
\references{
     M. Hubert, P. J. Rousseeuw, K. Vanden Branden (2005), ROBPCA: a new 
    approach to robust principal components analysis, \emph{Technometrics}, \bold{47}, 64--79.   
}
\seealso{
  \code{\link{FastHCS}}
}
\examples{
data(Tablets)
alpha<-0.5
Q<-15
p<-ncol(Tablets[,-1])
ns<-FHCSnumStarts(q=Q,eps=(1-alpha)*4/5)
RunExample<-FALSE
if(RunExample){
	Fit<-FastHCS(x=Tablets[,-1],q=Q,nSamp=ns,seed=1,alpha=0.5)
	colvec<-rep("orange",nrow(Tablets))
	colvec[Tablets[,1]==1]<-"blue"
	plot(Fit,col=colvec,pch=16)
}
}
\keyword{plot}
\keyword{robust}
\keyword{multivariate}
\author{Kaveh Vakili}
