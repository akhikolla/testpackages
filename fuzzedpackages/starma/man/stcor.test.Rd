\name{stcor.test}
\alias{stcor.test}
\alias{stcor.test.default}
\alias{print.stcor.test}
\title{Space-time series non correlation test}
\description{
\code{stcor.test} computes an extension of the Box-Pierce test statistic to accept or reject the
non correlation of the distinct observations of a given space-time series.
It is one of the main functions for the diagnostic part of the three-stage iterative model building procedure.
}
\usage{
stcor.test(data, wlist, tlag=NULL, slag=NULL, fitdf=0)

\method{print}{stcor.test}(x, ...)
}
\arguments{
\item{data}{ a matrix or data frame containing the space-time series: row-wise should be the temporal observations, with each column corresponding to a site. }
\item{wlist}{ a list of the weight matrices for each k-th order neighbours, first one being the identity. }
\item{tlag}{ the maximum time lag for the space-time autocorrelation functions. If \code{tlag = NULL}, it will use a large enough number of time lags. }
\item{slag}{ the maximum space lag for the space-time autocorrelation functions. If \code{slag = NULL}, it will use as many space lags as possible (as many as \code{length(wlist)}). }
\item{fitdf}{ number of degrees of freedom to be subtracted if \code{data} is a series of residuals. }
\item{x}{ a \code{starma} class object. }
\item{...}{ unused }
}
\details{
Since (Pfeifer and Deutsch, 1981) gives:
\deqn{Var(\hat{\rho}_l(s)) \approx \frac{1}{N(T-s)}}{Var(rho[l](s)) ~ 1 / N(T-s)}
We can extend Box-Pierce test statistic to space-time series:
\deqn{N \sum (T-s) \hat{\rho}_l(s)^2 \sim \chi^2(slag \times tlag)}{N \sum[(T-s)*rho[l](s)^2] ~ chisq(slag*tlag)}

\code{stcor.test} can be applied to a space-time series to test the null hypothesis of non correlation.

It is useful to check if the residuals of a STARMA models are multivariate white noise.
In this case, \code{fitdf} should be set equal to the number of parameters in the model.

Please note that this is an empirical extension and it has not yet been the subject of a paper.
The specifications of the weight matrices has not been studied either and could lead to inconsistencies.
}
\value{
A \code{data.frame} containing the following elements:
\item{X_squared}{ The value of the chi squared statistic }
\item{df}{ The degrees of freedom of the statistic (taking \code{fitdf} into account) }
\item{p.value}{ The p-value of the test }
}
\author{Felix Cheysson}
\references{
- Pfeifer, P., & Deutsch, S. (1980). A Three-Stage Iterative Procedure for Space-Time Modeling. Technometrics, 22(1): 35-47. doi:10.1080/00401706.1980.10486099

- Pfeifer, P., & Deutsch, S. (1981). Variance of the Sample Space-Time Autocorrelation Function. Journal of the Royal Statistical Society. Series B (Methodological), 43(1): 28-33.
}
\examples{
data(nb_mat)

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

# Test for multivariate normality
stcor.test(sim, blist)	# Data is correlated
stcor.test(eps, blist)		# Data should not be correlated (unless you're 5% unlucky)
}
\keyword{starma}
\keyword{estimation}
