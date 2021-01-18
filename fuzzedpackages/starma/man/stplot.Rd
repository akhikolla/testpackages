\name{stplot}
\alias{stplot}
\title{Plot for space-time series autocorrelation functions}
\description{
\code{stplot} renders a nice 2d plot for autocorrelation functions.
}
\usage{
stplot(acf, ci, call, ggplot=T)
}
\arguments{
\item{acf}{ a matrix containing the autocorrelation functions of a given space-time series: row-wise should be the temporal observations, with each column corresponding to a space lag. }
\item{ci}{ confidence intervals for the autocorrelation functions. }
\item{call}{ the name of the plot. }
\item{ggplot}{ a boolean indicating whether to use ggplot2 functions (they are recommended). }
}
\details{
This function plots the calculated autocorrelation functions of a space-time series.

In practice, the user should not use this function, as it is being called automatically when using \code{stacf} or \code{stpacf}.
The confidence intervals for the autocorrelation functions are approximated by 
\deqn{Var \left(\hat{\rho}_l(k)\right) \approx \frac{1}{N(T-k)}}{var(rho[lk]) = 1 / (N*(T-k))}
where N is the number of sites, and T the number of temporal observations.
}
\value{
NULL
}
\author{Felix Cheysson}
\references{
- Pfeifer, P., & Deutsch, S. (1981). Variance of the Sample Space-Time Autocorrelation Function. Journal of the Royal Statistical Society. Series B (Methodological), 43(1): 28-33.
}
\examples{
data(nb_mat)	# Get neighbourhood matrices

# Simulate a STARMA model
eps <- matrix(rnorm(94*200), 200, 94)
sim <- eps
for (t in 3:200) {
	sim[t,] <- (.4*diag(94) + .25*blist[[2]]) \%*\% sim[t-1,] +
		(.25*diag(94)                ) \%*\% sim[t-2,] +
		(            - .3*blist[[2]]) \%*\% eps[t-1,] +
		eps[t, ]
}

sim <- sim[101:200,]
sim <- stcenter(sim)	# Center and scale the dataset

# Autocorrelation functions
sim.stacf <- stacf(sim, blist, plot=FALSE)

# Plot the autocorrelation function
stplot(sim.stacf, 2 / sqrt(nrow(sim) * ncol(sim)), "stacf(sim, blist)")
}
\keyword{starma}
\keyword{plot}
