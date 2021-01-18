\name{FastHCS}
\alias{FastHCS}
\title{Performs the FastHCS algorithm for robust PCA.}
\usage{
  FastHCS(x,nSamp=NULL,alpha=0.5,q=10,seed=1)
}
\arguments{
  \item{x}{A numeric n (n>5*q) by p (p>1) matrix or data frame.}
   \item{nSamp}{A positive integer giving the number of resamples required;
          \code{"nsamp"} may not be reached if too many of the q-subsamples,
          chosen out of the observed vectors, are in a hyperplane.  If \code{"nSamp"} is
          omitted, it is calculated to provide a breakdown point of
          \code{"alpha"} with probability 0.99.}
 \item{alpha}{Numeric parameter controlling the size of the active subsets i.e., \code{"h=quanf(alpha,n,q)"}. Allowed
          values are between 0.5 and 1 and the default is 0.5.}
\item{q}{Number of principal components to compute. Note that p>q>1, 1<q<n. Default is q=10.}
\item{seed}{Starting value for random generator. Default is seed = 1.}
}
\value{
 A  list with components:
		\item{rawBest:}{The indexes of the h members of H*, the raw FastHCS optimal subset.}
	 	\item{obj:}{The FastHCS objective function corresponding to H*, the selected subset of h observations.}
	 	\item{rawDist:}{Outlyingness index of the data on the raw q-dimensonal subset that initialized H*.}
		\item{best:}{the indexes of the members of the H+, the FastHSC subset after the C-steps.} 
		\item{center:}{the p-vector of column means of the observations with indexes in \code{best}.}
		\item{loadings:}{the (rank q) loadings matrix of the observations with indexes in \code{best}.}
		\item{eigenvalues:}{the first \code{q)} eigenvalues of the observations with indexes in \code{best}.}
		\item{od:}{the orthogonal distances of the centered data wrt to the subspace spanned 
		by the \code{loadings} matrix.}
		\item{sd:}{the score distances of the data projected on the subspace spanned by the \code{loadings} matrix with respect to the estimated \code{center}.}	 
		\item{cutoff.od:}{the cutoff for the vector of orthogonal distances.}
		\item{cutoff.sd:}{the cutoff for the vector of score distances.}
		\item{scores}{The value of the projected on the space of the principal components data (the centred  data multiplied by the loadings matrix) is returned. Hence, cov(scores) is the diagonal matrix diag(eigenvalues).}
}
\description{
   Computes a robust PCA model with q components for 
	an n by p matrix of multivariate data using the FastHCS algorithm.
}
\examples{
## testing outlier detection
n<-100
p<-30
Q<-5
set.seed(123)
x0<-matrix(rnorm(n*p),nc=p)
x0[1:30,]<-matrix(rnorm(30*p,4.5,1/100),nc=p)
z<-c(rep(0,30),rep(1,70))
nStarts<-FHCSnumStarts(q=Q,eps=0.4)
Fit<-FastHCS(x=x0,nSamp=nStarts,q=Q)
z[Fit$best]
plot(Fit,col=(!z)+1,pch=16)

## testing outlier detection, different value of alpha
n<-100
p<-30
Q<-5
set.seed(123)
x0<-matrix(rnorm(n*p),nc=p)
x0[1:20,]<-matrix(rnorm(20*p,4.5,1/100),nc=p)
z<-c(rep(0,20),rep(1,80))
nStarts<-FHCSnumStarts(q=Q,eps=0.25)
Fit<-FastHCS(x=x0,nSamp=nStarts,q=Q,alpha=0.75)
z[Fit$best]


#testing exact fit
n<-100
p<-5
Q<-4
set.seed(123)
x0<-matrix(rnorm(n*p),nc=p)
x0[1:30,]<-matrix(rnorm(30*p,4.5,1/100),nc=p)
x0[31:100,4:5]<-x0[31:100,2]
z<-c(rep(0,30),rep(1,70))
nStart<-FHCSnumStarts(q=Q,eps=0.4)
results<-FastHCS(x=x0,nSamp=nStart,q=Q)
z[results$best]
results$obj

#testing rotation equivariance
n<-100
p<-10
Q<-3
set.seed(123)
x0<-scale(matrix(rnorm(n*p),nc=p))
A<-diag(rep(1,p))
A[1:2,1:2]<-c(0,1,-1,0)
x1<-x0\%*\%A
nStart<-FHCSnumStarts(q=Q,eps=0.4)
r0<-FastHCS(x=x0,nSamp=nStart,q=Q,seed=0)
r1<-FastHCS(x=x1,nSamp=nStart,q=Q,seed=0)
max(abs(log(r1$eigenvalues[1:Q]/r0$eigenvalues[1:Q])))
}
\references{
	 Schmitt E. and Vakili K. and (2015).
	 Robust PCA with FastHCS.
	(http://arxiv.org/abs/1402.3514)
}
\author{
	Kaveh Vakili, Eric Schmitt
}
\keyword{multivariate}
\keyword{robust}

