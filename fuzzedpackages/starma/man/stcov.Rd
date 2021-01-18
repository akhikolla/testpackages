\name{stcov}
\alias{stcov}
\title{Space-time covariance function}
\description{
\code{stcov} computes the space-time covariance of the serie \code{data} between \code{slag1}-th and \code{slag2}-th order neighbours at time lag \code{tlag}.
}
\usage{
stcov(data, wlist, slag1, slag2, tlag)
}
\arguments{
\item{data}{ a matrix or data frame containing the space-time series: row-wise should be the temporal observations, with each column corresponding to a site. }
\item{wlist}{ a list of the weight matrices for each k-th order neighbours, first one being the identity. }
\item{slag1, slag2}{ the space lags for the space-time covariance. }
\item{tlag}{ the time lag for the space-time covariance. }
}
\details{
\code{stcov} is mainly used as an internal function for the computation of \code{\link{stacf}} and \code{\link{stpacf}}.
\code{slag1} and \code{slag2} must be lower than \code{length(wlist)}.

It is computed as follows:
\deqn{\hat{\gamma}_{lk}(s) = \frac{1}{N(T-s)} Tr \left( \sum_{t=s+1}^{T} W^{(k)\prime} W^{(l)} z_t z_{t-k}^\prime \right)}{gamma[lk](s) = Tr ( \sum_[t=s+1]^[T] ( t(W[k]) * W[l] * z(t) * t(z(t-k)) )) / N(T-s) }
}
\value{
A numeric.
}
\author{Felix Cheysson}
\references{
Pfeifer, P., & Deutsch, S. (1980). A Three-Stage Iterative Procedure for Space-Time Modeling. Technometrics, 22(1), 35-47. doi:10.1080/00401706.1980.10486099
}
\examples{
data(nb_mat)	# Get neighbourhood matrices

data <- matrix(rnorm(9400), 100, 94)

# Compute covariance between 2-nd and 1-st order neighbours, at time lag 5
stcov(data, blist, 2, 1, 5)		
}
\keyword{starma}
\keyword{stcov}
