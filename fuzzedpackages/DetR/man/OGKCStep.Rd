\name{OGKCStep}
\alias{OGKCStep}
\title{Robust and Deterministic Linear Regression via OGKCStep}
\description{
     Function to find the OGKCStep ('best') H-subset.
}
\usage{
  OGKCStep(x0, scale_est, alpha=0.5)
}
\arguments{
  \item{x0}{Matrix of continuous variables. }
  \item{alpha}{numeric parameter controlling the size of the subsets over which the determinant is minimized, i.e., alpha*n observations are used for computing the determinant. Allowed values are between 0.5 and 1 and the default is 0.5.}
  \item{scale_est}{A character string specifying the
  variance functional.  Possible values are Qn or scaleTau2.}
}
\value{
 
  \item{best}{
    the best subset found and used for computing the raw estimates, with
    \code{\link{length}(best) == quan = \link{h.alpha.n}(alpha,n,p)}.
  }
 

}
\references{
    Maronna, R.A. and Zamar, R.H. (2002) Robust estimates of location
     and dispersion of high-dimensional datasets; \emph{Technometrics}
     \bold{44}(4), 307--317.

 Rousseeuw, P.J. and Croux, C. (1993) Alternatives to the Median Absolute Deviation; 
  \emph{Journal of the American Statistical Association} , \bold{88}(424), 1273--1283.

  Peter J. Rousseeuw (1984), Least Median of Squares Regression.
  \emph{Journal of the American Statistical Association} \bold{79}, 871--881.

  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection.} Wiley.

  P. J. Rousseeuw and K. van Driessen (1999)
  A fast algorithm for the minimum covariance determinant estimator.
  \emph{Technometrics} \bold{41}, 212--223.

  Pison, G., Van Aelst, S., and Willems, G. (2002)
  Small Sample Corrections for LTS and MCD.
  \emph{Metrika} \bold{55}, 111--123.

  Hubert, M., Rousseeuw, P. J. and Verdonck, T. (2012) A deterministic algorithm for robust location
 and scatter. \emph{Journal of Computational and Graphical Statistics}
  \bold{21}, 618--637.
}
\examples{
n<-100
set.seed(123)# for reproducibility
x0<-matrix(rnorm(n*2),nc=2)
out1<-OGKCStep(x0,alpha=0.5,scale_est=pcaPP::qn)

#comparaison with DetMCD:

#a) create data

set.seed(123456)
Simulation<-DetR:::fx01()
#should be \approx 10
sqrt(min(mahalanobis(Simulation$Data[Simulation$label==0,],rep(0,ncol(Simulation$Data)),
Simulation$Sigma_u))/qchisq(0.975,df=ncol(Simulation$Data)))
a0<-eigen(Simulation$Sigma_u)
Su_ih<-(a0$vector)\%*\%diag(1/sqrt(a0$values))\%*\%t(a0$vector)
#run algorithms 
A0<-robustbase::covMcd(Simulation$Data,nsamp='deterministic',scalefn=pcaPP::qn,alpha=0.5)
A1<-OGKCStep(Simulation$Data,alpha=0.5,scale_est=pcaPP::qn)
#getbiases algorithms 
SB<-eigen(Su_ih\%*\%var(Simulation$Data[A1,])\%*\%Su_ih)$values
log10(SB[1]/SB[ncol(Simulation$Data)-1])
SB<-eigen(Su_ih\%*\%var(Simulation$Data[A0$best,])\%*\%Su_ih)$values
log10(SB[1]/SB[ncol(Simulation$Data)-1])
}
\author{
  Large part of the the code are from function \code{.detmcd} in package robustbase , 
 , see citation("robustbase")
}
\keyword{multivariate}
\keyword{robust}
\keyword{deterministic}
