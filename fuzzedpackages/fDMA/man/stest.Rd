\name{stest}
\alias{stest}

\title{Computes a Few Stationarity Tests.}

\description{
This is a wrapper for three functions from \code{tseries} package. Augmented Dickey-Fuller (ADF, \code{\link[tseries]{adf.test}}), Phillips-Perron (PP, \code{\link[tseries]{pp.test}}) and Kwiatkowski-Phillips-Schmidt-Shin (KPSS, \code{\link[tseries]{kpss.test}}) tests for stationarity are performed. 
}

\usage{
stest(data)
}

\arguments{
\item{data}{\code{\link[base]{matrix}} of variables, different columns correspond to different variables}
}

\value{\code{\link[base]{matrix}}, 
tests statistics and p-values are given by columns, 
tests outcomes for different variables are ordered by rows}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
ld.wti <- (diff(log(wti)))[-1,]
ld.drivers <- (diff(log(drivers)))[-1,]
x <- cbind(ld.wti,ld.drivers)
stest(x)
}
}
