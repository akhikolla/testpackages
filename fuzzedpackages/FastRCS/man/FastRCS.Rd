\name{FastRCS}
\alias{FastRCS}
\title{Computes the FastRCS outlyingness index for regression.}
\usage{
  FastRCS(x,y,nSamp,alpha=0.5,seed=1,intercept=1)
}
\arguments{
  \item{x}{A numeric n (n>5*p) by p (p>1) matrix or data frame. Should 
	not contain an intercept.}
  \item{y}{A numeric nvector.}
   \item{nSamp}{a positive integer giving the number of resamples required;
          \code{"nSamp"} may not be reached if too many of the p-subsamples,
          chosen out of the observed vectors, are in a hyperplane.  If \code{"nSamp"} is
          omitted, it is calculated so that the probability of getting at least one
      uncontaminated starting point is always at least 99 percent when there are n/2 outliers.}
 \item{alpha}{numeric parameter controlling the size of the active subsets,
	 i.e., \code{"h=quanf(alpha,n,p)"}.  Allowed
          values are between 0.5 and 1 and the default is 0.5.}
\item{seed}{starting value for random generator. A positive integer. Default is seed = 1}
	\item{intercept}{If true, a model with constant term will be estimated; otherwise no constant term will be included. Default is intercept=TRUE.}

}
\value{
\item{nSamp}{The value of nSamp used.}
\item{alpha}{The value of alpha used.}
\item{obj}{The value of the FastRCS objective function (the I-index) obtained for H*.}
\item{rawBest}{The index of the h observation with smallest outlyingness indexes.}
\item{rawDist}{The distances of the observations to the model defined by rawBest.}
\item{best}{The index of the J observation with outlyingness smaller than the rejection threshold.}
\item{coefficients}{The vector of coefficients of the hyperplane fitted to the members of \code{$rew$best}.}
\item{fitted.values}{the fitted mean values: \code{cbind(1,x)\%*\%rew$coefficients}.}
\item{residuals}{the residuals, that is response minus fitted values.}
\item{rank}{the numeric rank of the fitted linear model.}
\item{weights}{(only for weighted fits) the specified weights.}
\item{df.residual}{the residual degrees of freedom.}
\item{scale}{(robust) scale estimate of the reweighted residuals.}
}

\details{
	The current version of FastRCS includes the use of a C-step procedure to improve efficiency (Rousseeuw and van Driessen (1999)). C-steps are taken after the raw subset is found and before reweighting. In experiments, we found that carrying C-Steps 
	starting from the members of \code{$rawBest} improves the speed of convergence without increasing the bias
	 of the final estimates. FastRCS is regression and affine equivariant and thus consistent at the elliptical
	 model (Grubel and Rock (1990)). 
}
\description{
   Computes a fast and robust regression model for a n by p matrix of 
	multivariate continuous regressors and a single dependent variable.
}
\examples{
## testing outlier detection
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
z[results$best]

## testing outlier detection, different value of alpha
set.seed(123)
n<-100
p<-3
x0<-matrix(rnorm(n*p),nc=p)
y0<-rnorm(n)
z<-c(rep(0,20),rep(1,80))
x0[1:20,]<-matrix(rnorm(20*p,5,1/100),nc=p)
y0[1:20]<-rnorm(20,10,1/100)
ns<-FRCSnumStarts(p=p,eps=0.25);
results<-FastRCS(x=x0,y=y0,alpha=0.75,nSamp=ns)
z[results$best]



#testing exact fit
set.seed(123)
n<-100
p<-3
x0<-matrix(rnorm(n*p),nc=p)
y0<-rep(1,n)
z<-c(rep(0,30),rep(1,70))
x0[1:30,]<-matrix(rnorm(30*p,5,1/100),nc=p)
y0[1:30]<-rnorm(30,10,1/100)
ns<-FRCSnumStarts(p=p,eps=0.4);
results<-FastRCS(x=x0,y=y0,alpha=0.5,nSamp=ns,seed=1)
z[results$rawBest]
results$obj

#testing regression equivariance
n<-100
p<-3
x0<-matrix(rnorm(n*(p-1)),nc=p-1)
y0<-rnorm(n)
ns<-FRCSnumStarts(p=p,eps=0.4);
y1<-y0+cbind(1,x0)\%*\%rep(-1,p)
results1<-FastRCS(y=y0,x=x0,nSamp=ns,seed=1)$coefficients
results2<-FastRCS(y=y1,x=x0,nSamp=ns,seed=1)$coefficients
results1+rep(-1,p)
#should be the same:
results2
}
\references{
	Grubel, R. and Rocke, D. M. (1990).
	On the cumulants of affine equivariant estimators in elliptical families.
	Journal of Multivariate Analysis, Vol. 35, p. 203--222. Journal of Multivariate Analysis

	Rousseeuw, P. J., and van Driessen, K. (2006). Computing lts regression for large data sets. 
	Data mining and	Knowledge Discovery, 12, 29--45

	Vakili, K. and Schmitt, E. (2014).
	Finding Regression Outliers With FastRCS.
	(http://arxiv.org/abs/1307.4834)
}
\author{
	Kaveh Vakili
}
\keyword{regression}
\keyword{robust}

