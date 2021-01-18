\name{MultipleFeatures}
\alias{MultipleFeatures}
\title{Fourier coefficients describing the shape of many hand written replications of the numerals '0' and '1'.}
\description{
  A data frame with the subset of the 'Multiple Features'
 dataset corresponding to the sample of '0' and '1' numerals. 
}
\usage{
data(MultipleFeatures)
}
\format{
  \describe{
	\item{Labels}{Numerals.}
	\item{Column 2--77}{Fourier coefficients describing the shape of each observation.}
  }
}
\source{
Van Breukelen, M. Duin, R.P.W. Tax, D.M.J. and  Den Hartog, J.E. (1998).
Handwritten digit recognition by combined classifiers, Kybernetika, vol. 34, 381--386.
}
\examples{
data(MultipleFeatures)
alpha<-0.5
Q<-15
p<-ncol(MultipleFeatures[,-1])
ns<-FHCSnumStarts(q=Q,eps=(1-alpha)*4/5)
RunExample<-FALSE
if(RunExample){
	Fit<-FastHCS(x=MultipleFeatures[,-1],q=Q,nSamp=ns,seed=1)
	colvec<-rep("orange",nrow(MultipleFeatures))
	colvec[MultipleFeatures[,1]==1]<-"blue"
	plot(Fit,col=colvec,pch=16)
}

}
\keyword{datasets}
