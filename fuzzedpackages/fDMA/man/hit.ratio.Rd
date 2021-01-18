\name{hit.ratio}
\alias{hit.ratio}

\title{Computes Hit Ratio (HR) for Forecast.}

\description{
Sometimes it is interesting to analyze just whether the forecast can predict the direction of a change in a modelled time-series. This function computes the proportion of correctly predicted signs (i.e., in which cases the direction of a change given by forecast agrees with the change in real data). 
}
\usage{
hit.ratio(y,y.hat,d=NULL)
}

\arguments{
\item{y}{\code{\link[base]{numeric}}, \code{\link[base]{vector}}, or one row or one column \code{\link[base]{matrix}} or \code{\link[xts]{xts}} object, representing a forecasted time-series}
\item{y.hat}{\code{\link[base]{numeric}}, \code{\link[base]{vector}}, or one row or one column \code{\link[base]{matrix}} or \code{\link[xts]{xts}} object, representing forecast predictions}
\item{d}{optional, \code{\link[base]{logical}}, \code{d=FALSE} for level time-series, \code{d=TRUE} if time-series already represent changes, by default \code{d=FALSE}}
}

\value{
\code{\link[base]{numeric}}
}

\examples{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
\dontrun{
m1 <- fDMA(y=wti,x=drivers,alpha=0.99,lambda=0.99,initvar=10)
hit.ratio(y=wti,y.hat=m1$y.hat)

m2 <- fDMA(y=ld.wti,x=ld.drivers,alpha=0.99,lambda=0.99,initvar=10)
hit.ratio(y=ld.wti,y.hat=m2$y.hat,d=TRUE)
}
\dontshow{
m1 <- fDMA(y=wti[1:20,1],x=drivers[1:20,1:4],alpha=0.99,lambda=0.99,initvar=10)
hit.ratio(y=wti[1:20,1],y.hat=m1$y.hat)

m2 <- fDMA(y=ld.wti[1:20,1],x=ld.drivers[1:20,1:4],alpha=0.99,lambda=0.99,initvar=10)
hit.ratio(y=ld.wti[1:20,1],y.hat=m2$y.hat,d=TRUE)
}
}

\references{
Baur, D. G., Beckmann, J., Czudaj, R., 2016. A melting pot -- Gold price forecasts under model and parameter uncertainty. \emph{International Review of Financial Analysis} \bold{48}, 282--291. 
}
