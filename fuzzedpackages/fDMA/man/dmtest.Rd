\name{dmtest}
\alias{dmtest}

\title{Computes Diebold-Mariano Test.}

\description{
This is a wrapper for \code{\link[forecast]{dm.test}} from \code{forecast} package. This function computes the original Diebold-Mariano test. 
}

\details{The null hypothesis is that the two methods have the same forecast accuracy. This function assumes that one-step ahead forecasts are compared and the second power is used in the loss function (see \code{\link[forecast]{dm.test}}). "The Diebold-Mariano (DM) test was intended for comparing forecasts; it has been, and remains, useful in that regard. The DM test was not intended for comparing models." (Diebold, 2015)}

\usage{
dmtest(y,f)
}

\arguments{
\item{y}{\code{\link[base]{vector}} of the forecasted time-series}
\item{f}{\code{\link[base]{matrix}} of the predicted values from various methods, forecasts are ordered in rows, the first row should correspond to the method that is compared with alternative ones (corresponding to subsequent rows)}
}

\value{\code{\link[base]{matrix}}, 
first column contains tests statistics, next p-values are given for the alternative hypothesis that alternative forecasts have different accuracy than the compared forecast, alternative forecasts are less accurate and alternative forecasts have greater accuracy, tests outcomes for different forecasts are ordered by rows}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
m <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.90,initvar=10)
m <- m$y.hat
a <- altf2(y=ld.wti,x=ld.drivers,d=TRUE)
a <- a$y.hat
a <- matrix(unlist(a),nrow=length(a),byrow=TRUE)
fc <- rbind(m,a)
dm <- dmtest(y=as.vector(ld.wti),f=fc)
}
}

\references{
Diebold, F. X., 2015. Comparing predictive accuracy, Twenty years later: A peersonal perspective on the use and abuse of Diebold-Mariano tests. \emph{Journal of Business & Economic Statistics} \bold{33}, \url{http://dx.doi.org/10.1080/07350015.2014.983236}.

Diebold, F. X., Mariano, R. S., 1995. Comparing predictive accuracy. \emph{Journal of Business & Economic Statistics} \bold{13}, 253--263.
}

\seealso{\code{\link{hmdmtest}}, \code{\link{mdmtest}}.}
