\name{tvp}
\alias{tvp}

\title{Computes Time-Varying Parameters Regression.}

\description{
This function computes Time-Varying Parameters Regression (TVP) with the updating procedure as in Raftery et. al (2010). 
}

\usage{
tvp(y,x,V,lambda,W=NULL,kappa=NULL,c=NULL)
}

\arguments{
\item{y}{\code{\link[base]{numeric}} or a column \code{\link[base]{matrix}} of a dependent variable}
\item{x}{\code{\link[base]{matrix}} of independent variables, different columns should correspond to different variables}
\item{V}{\code{\link[base]{numeric}}, initial variance in the state space equation for the recursive moment estimator updating method, as in Raftery et al. (2010)}
\item{lambda}{\code{\link[base]{numeric}}, a forgetting factor between 0 and 1 used in variance approximations}
\item{W}{optional, \code{\link[base]{numeric}}, initial value of variance for the model equations, if not specified the method based on the linear regression, as in Raftery et al. (2010) is used}
\item{kappa}{optional, \code{\link[base]{numeric}}, a parameter in the exponentially weighted moving average in variance updating (see also \code{\link{fDMA}}), between 0 and 1, if not specified the method as in Raftery et al. (2010) is used}
\item{c}{optional, \code{\link[base]{logical}}, a parameter indicating whether constant is included, if not specified \code{c=TRUE} is used, i.e., constant is included}
}

\details{It is not possible to set \code{c=FALSE} if \code{ncol(x)=0}. In such a case the function will automatically reset \code{c=TRUE} inside the code.}

\value{class \code{tvp} object, \code{\link[base]{list}} of

\item{$y.hat}{fitted (forecasted) values}
\item{$thetas}{estimated regression coefficients}
\item{$pred.dens.}{predicitive densities from each period}
\item{$y}{\code{y}, forecasted time-series}
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
t1 <- tvp(y=ld.wti,x=ld.drivers,V=1,lambda=0.99)
}
\dontshow{
tt1 <- tvp(y=ld.wti[1:20,1],x=ld.drivers[1:20,1:4],V=1,lambda=0.99)
}
\dontrun{
t2 <- tvp(y=ld.wti,x=ld.drivers,V=1,lambda=0.99,W=1)
}
\dontshow{
tt2 <- tvp(y=ld.wti[1:20,1],x=ld.drivers[1:20,1:4],V=1,lambda=0.99,W=1)
}
\dontrun{
t3 <- tvp(y=ld.wti,x=ld.drivers,V=1,lambda=0.99,W=1,kappa=0.75)
}
\dontshow{
tt3 <- tvp(y=ld.wti[1:20,1],x=ld.drivers[1:20,1:4],V=1,lambda=0.99,W=1,kappa=0.75)
}
# Model with constant only
empty <- matrix(,nrow=nrow(ld.drivers),ncol=0)
t4 <- tvp(y=ld.wti,x=empty,lambda=0.99,V=1)
}

\references{
Raftery, A. E., Karny, M., Ettler, P., 2010. Online prediction under model uncertainty via Dynamic Model Averaging: Application to a cold rolling mill. \emph{Technometrics} \bold{52}, 52--66.

Sanderson, C., Curtin, R., 2016. Armadillo: A template-based C++ library for linear algebra. \emph{Journal of Open Source Software} \bold{1}, \url{http://arma.sourceforge.net/armadillo_joss_2016.pdf}.
}

\seealso{\code{\link{grid.tvp}}, \code{\link{print.tvp}}, \code{\link{summary.tvp}}, \code{\link{plot.tvp}}.}
