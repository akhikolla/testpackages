\name{archtest}
\alias{archtest}

\title{Computes Engle's ARCH Test.}

\description{
This function computes Engle's ARCH test. The null hypothesis of this Lagrange Multiplier test is that a series of residuals exhibits no ARCH effects. The alternative hypothesis is that ARCH(lag) effects are present. The \code{lag} is specified by the User.
}

\usage{
archtest(ts,lag=NULL)
}

\arguments{
\item{ts}{\code{\link[base]{vector}}, the tested time-series}
\item{lag}{\code{\link[base]{numeric}}, suspected order of ARCH process, if not specified \code{lag=1} is taken}
}

\value{class \code{htest} object, \code{\link[base]{list}} of 
\item{statistic}{test statistic}
\item{parameter}{\code{lag} used in the test}
\item{alternative}{alternative hypothesis of the test}
\item{p.value}{p-value}
\item{method}{name of the test}
\item{data.name}{name of the tested time-series}
}

\examples{
wti <- crudeoil[-1,1]
ld.wti <- (diff(log(wti)))[-1,]
arch <- archtest(ts=as.vector(ld.wti),lag=10)
}

\references{
Engle, R. F., 1982. Autoregressive conditional heteroscedasticity with estimates of the variance of United Kingdom inflation. \emph{Econometrica} \bold{50}, 987--1007.
}
