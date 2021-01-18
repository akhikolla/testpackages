\name{tvp.reg}
\alias{tvp.reg}

\title{Computes Time-Varying Parameters Regression.}

\description{
This function estimates Time-Varying Parameters regression. 
}

\usage{
tvp.reg(y,x,lambda=NULL,kappa=NULL,V=NULL,W=NULL)
}

\arguments{
\item{y}{one column {\code{\link[base]{matrix}}} of forecasted time-series, observations inserted rowwise}
\item{x}{{\code{\link[base]{matrix}}} of independent time-series (predictors), observations inserted rowwise}
\item{lambda}{optional, see \code{\link[dynmix]{mixest1}}}
\item{kappa}{optional, see \code{\link[dynmix]{mixest1}}}
\item{V}{optional, \code{\link[base]{numeric}} initial variance, by default \code{V=1} is taken}
\item{W}{optional, \code{\link[base]{numeric}} initial value to be put on diagonal of covariance matrix, by default \code{W=1} is taken}
}

\details{
If \code{lambda} is specified, then the method described by Raftery et al. (2010) is used, with possible extentsion to the one described by Koop and Korobilis (2012). Otherwise, the Kalman filter described as by Nagy and Suzdaleva (2013) is used.  
}

\value{object of class {\code{tvpreg}}, i.e., {\code{\link[base]{list}}} of 
\item{$y.hat}{{\code{\link[base]{vector}}} of predictions}
\item{$coef}{{\code{\link[base]{matrix}}} of regression coefficients}
\item{$R}{{\code{\link[base]{matrix}}} of diagonals of covariances corresponding to independent variables in regressions}
\item{$V}{{\code{\link[base]{vector}}} of outcome variances}
}

\examples{
data(oil)
t1 <- tvp.reg(y=oil[,1,drop=FALSE],x=oil[,-1,drop=FALSE],lambda=0.99,V=100,W=100)
t2 <- tvp.reg(y=oil[,1,drop=FALSE],x=oil[,-1,drop=FALSE],V=100,W=100)
}

\references{
Koop, G., Korobilis, D., 2012, Forecasting inflation using Dynamic Model Averaging. \emph{International Economic Review} \bold{53}, 867--886.

Nagy, I., Suzdaleva, E., 2017, \emph{Algorithms and Programs of Dynamic Mixture Estimation}, Springer. 

Raftery, A. E., Karny, M., Ettler, P., 2010, Online prediction under model uncertainty via Dynamic Model Averaging: Application to a cold rolling mill. \emph{Technometrics} \bold{52}, 52--66.
}
