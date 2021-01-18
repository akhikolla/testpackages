\name{FastPCS}
\alias{FastPCS}
\title{Computes the FastPCS multivariate outlyingness index.}
\usage{
  FastPCS(x,nSamp,alpha=0.5,seed=1)
}
\arguments{
  \item{x}{A numeric n (n>5*p) by p (p>1) matrix or data frame.}
 \item{nSamp}{A positive integer giving the number of resamples required;
          \code{"nSamp"} may not be reached if too many of the p-subsets,
          chosen out of the observed vectors, lie on a hyperplane.  If \code{"nSamp"} is
          omitted, it is calculated so that the probability of getting at least one
      uncontaminated starting point is always at least 99 percent when there are n/2 outliers.}
 \item{alpha}{Numeric parameter controlling the size of the active subsets,
	 i.e., \code{"h=quanf(alpha,n,p)"}.  Allowed
          values are between 0.5 and 1 and the default is 0.5.}
\item{seed}{Starting value for random generator. A positive integer. Default is seed = 1}
}
\value{
\item{alpha}{The value of alpha used.}
\item{nSamp}{The value of nSamp used.}
\item{obj}{The value of the FastPCS objective function of the optimal h subset.}
\item{rawBest}{The index of the h observation with smallest outlyingness indexes.}
item{best}{The index of the observations with outlyingness smaller than the rejection threshold after C-steps are taken.}
\item{center}{The mean vector of the observations with outlyingness smaller than the rejection threshold after C-steps are taken.}
\item{cov}{Covariance matrix of the observations with outlyingness smaller than the rejection threshold after C-steps are taken.}
\item{distance}{The statistical distance of each observation wrt the center vector and cov matrix of the observations with outlyingness smaller than the rejection threshold after C-steps are taken.}
}
\details{

	The current version of FastPCS includes the use of a C-step procedure to improve efficiency (Rousseeuw and van Driessen (1999)). C-steps are taken after the raw subset (H*) as been chosen (according to the I-index) and before reweighting. In experiments, we found that carrying C-Steps 
	starting from the members of \code{$rawBest} improves the speed of convergence without increasing the bias
	 of the final estimates. FastPCS is affine equivariant (Schmitt et al. (2014)) and thus consistent at the 
	elliptical model (Maronna et al., (2006) p. 217).
}
\description{
   Computes a fast and robust multivariate outlyingness index for a n by p matrix of 
	multivariate continuous data.
}
\examples{
## testing outlier detection 
set.seed(123)
n<-100
p<-3
x0<-matrix(rnorm(n*p),nc=p)
x0[1:30,]<-matrix(rnorm(30*p,4.5,1/100),nc=p)
z<-c(rep(0,30),rep(1,70))
nstart<-FPCSnumStarts(p=p,eps=0.4)
results<-FastPCS(x=x0,nSamp=nstart)
z[results$best]

## testing outlier detection, different value of alpha
set.seed(123)
n<-100
p<-3
x0<-matrix(rnorm(n*p),nc=p)
x0[1:20,]<-matrix(rnorm(20*p,4.5,1/100),nc=p)
z<-c(rep(0,20),rep(1,80))
nstart<-FPCSnumStarts(p=p,eps=0.25)
results<-FastPCS(x=x0,nSamp=nstart,alpha=0.75)
z[results$best]

#testing exact fit
set.seed(123)
n<-100
p<-3
x0<-matrix(rnorm(n*p),nc=p)
x0[1:30,]<-matrix(rnorm(30*p,5,1/100),nc=p)
x0[31:100,3]<-x0[31:100,2]*2+1
z<-c(rep(0,30),rep(1,70))
nstart<-FPCSnumStarts(p=p,eps=0.4)
results<-FastPCS(x=x0,nSamp=nstart)
z[results$rawBest]
results$obj

#testing affine equivariance
n<-100
p<-3
set.seed(123)
x0<-matrix(rnorm(n*p),nc=p)
nstart<-500
results1<-FastPCS(x=x0,nSamp=nstart,seed=1)
a1<-matrix(0.9,p,p)
diag(a1)<-1
x1<-x0\%*\%a1
results2<-FastPCS(x=x1,nSamp=nstart,seed=1)
results2$center
results2$cov
#should be the same
results1$center\%*\%a1
a1%*%results1$cov\%*\%t(a1)

}
\references{
	Maronna, R. A., Martin R. D. and Yohai V. J. (2006). Robust Statistics: Theory and Methods. Wiley, New York.

	P. J. Rousseeuw and K. van Driessen (1999). A fast algorithm for the minimum covariance determinant estimator. Technometrics 41, 212--223.

	Eric Schmitt, Viktoria Oellerer, Kaveh Vakili (2014).
	The finite sample breakdown point of PCS Statistics and Probability Letters, Volume 94, Pages 214-220.

	Vakili, K. and Schmitt, E. (2014).
	Finding multivariate outliers with FastPCS.
	\emph{Computational Statistics \& Data Analysis}.
  	Vol. 69, pp 54--66.
	(http://arxiv.org/abs/1301.2053)
}
\author{
	Kaveh Vakili
}
\keyword{multivariate}
\keyword{robust}

