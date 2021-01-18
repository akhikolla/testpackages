\name{stacf}
\alias{stacf}
\alias{stpacf}
\title{Space-time autocorrelation functions}
\description{
The functions defined below are the main tools to the identification and the diagnostic part of the three-stage iterative model procedure building.
\code{stacf} and \code{stpacf} respectively compute the autocorrelation and partial autocorrelation functions of a space-time series.
}
\usage{
stacf(data, wlist, tlag.max=NULL, plot=TRUE, use.ggplot=TRUE)
stpacf(data, wlist, tlag.max=NULL, plot=TRUE, use.ggplot=TRUE)
}
\arguments{
\item{data}{ a matrix or data frame containing the space-time series: row-wise should be the temporal observations, with each column corresponding to a site. }
\item{wlist}{ a list of the weight matrices for each k-th order neighbours, first one being the identity. }
\item{tlag.max}{ the maximum time lag for the space-time autocorrelation functions. If \code{tlag.max = NULL}, it will use a large enough number of time lags. }
\item{plot}{ whether to plot the autocorrelation functions or not. }
\item{use.ggplot}{ if \code{plot = TRUE}, whether to use \code{ggplot2} or not to display the autocorrelation functions. Not using \code{ggplot2} is depreciated. }
}
\details{
\code{stacf} and \code{stpacf} respectively compute the space-time autocorrelation and partial autocorrelation functions of the serie \code{data} between \code{s}-th and \code{0}-th order neighbors at time lag \code{t}, 
for \code{s} ranging from \code{0} to \code{length(wlist)} and \code{t} ranging from \code{1} to \code{tlag.max}.

The autocorrelation function is computed as follows:
\deqn{\hat{\rho}_{l}(s) = frac{\hat{\gamma}_{l0}(s)}{\left( \hat{\gamma}_{ll}(0) \hat{\gamma}_{00}(s) \right) ^{1/2} }}{rho[l](s) = gamma[l0](s) / sqrt( (gamma[ll](0) * gamma[00](0)) )}

The partial autocorrelation functions are computed solving iteratively the Yule Walker equations for increasing time lags and space lags.

Note that the identification might be biased if the partial autocorrelation functions are not computed with enough space lags, since Yule Walker equations are sensible to the maximum space lag given.
}
\value{
An object of class \code{matrix} containing the estimated acf.
Row-wise are the different time lags, column-wise the different space lags.
}
\author{Felix Cheysson}
\references{
Pfeifer, P., & Deutsch, S. (1980). A Three-Stage Iterative Procedure for Space-Time Modeling. Technometrics, 22(1), 35-47. doi:10.1080/00401706.1980.10486099
}
\examples{
data(nb_mat)	# Get neighbourhood matrices

# Simulate a STARMA model
eps <- matrix(rnorm(94*200), 200, 94)
sim <- eps
for (t in 3:200) {
	sim[t,] <- (.4*blist[[1]] + .25*blist[[2]]) \%*\% sim[t-1,] +
		(.25*blist[[1]]                ) \%*\% sim[t-2,] +
		(            - .3*blist[[2]]) \%*\% eps[t-1,] +
		eps[t, ]
}

sim <- sim[101:200,]
sim <- stcenter(sim)	# Center and scale the dataset

# Plot stacf and stpacf
stacf(sim, blist)
stpacf(sim, blist)
}
\keyword{starma}
\keyword{identification}
\keyword{diagnostic}
\keyword{acf}
\keyword{stacf}
\keyword{pacf}
\keyword{stpacf}
