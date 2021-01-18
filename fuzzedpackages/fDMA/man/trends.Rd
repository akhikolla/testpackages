\name{trends}
\alias{trends}

\title{Google Trends for Crude Oil Data.}

\description{Google Trends for Crude Oil Data.}

\usage{data(trends)}

\format{\code{trends} is \code{\link[xts]{xts}} object such that
\itemize{
\item \code{trends$stock_markets} -- Google Trends for "stock markets"
\item \code{trends$interest_rate} -- Google Trends for "interest rate"
\item \code{trends$economic_activity} -- Google Trends for "economic activity"
\item \code{trends$exchange_rate} -- Google Trends for "exchange rate"
\item \code{trends$oil_production} -- Google Trends for "oil production"
\item \code{trends$oil_consumption} -- Google Trends for "oil consumption"
\item \code{trends$market_stress} -- Google Trends for "market stress"
}
}

\details{The data are in monthly frequency. They cover the period between Jan, 2004 and Dec, 2016.
}

\source{The data are provided by Google. 

\url{https://trends.google.com/trends}
}

\examples{
data(trends)
gtrends <- trends/100
data(crudeoil)
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
ld.wti <- ld.wti['2004-01-01-/']
ld.drivers <- ld.drivers['2004-01-01-/']
\dontrun{
xx <- ld.drivers
m <- fDMA(y=ld.wti,x=xx,alpha=0.99,lambda=0.99,initvar=1,model="dma",gprob=gtrends,omega=0.5)
}
}
