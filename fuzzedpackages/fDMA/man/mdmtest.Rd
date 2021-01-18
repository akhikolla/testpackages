\name{mdmtest}
\alias{mdmtest}

\title{Computes Harvey-Leybourne-Newbold Test.}

\description{
This is a wrapper for \code{\link[forecast]{dm.test}} from \code{forecast} package. This function computes the modified Diebold-Mariano test. The modification is useful for small samples. 
}

\details{The null hypothesis is that the two methods have the same forecast accuracy. This function assumes that one-step ahead forecasts are compared and the second power is used in the loss function (see \code{\link[forecast]{dm.test}}).}

\usage{
mdmtest(y,f)
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
mdm <- mdmtest(y=as.vector(ld.wti),f=fc)
}
}

\references{
Harvey, D., Leybourne, S., Newbold, P., 1997. Testing the equality of prediction mean squared errors. \emph{International Journal of Forecasting} \bold{13}, 281--291. 
}

\seealso{\code{\link{dmtest}}, \code{\link{hmdmtest}}.}
