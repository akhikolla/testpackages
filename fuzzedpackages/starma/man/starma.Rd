\name{starma}
\alias{starma}
\alias{starma.default}
\alias{print.starma}
\title{Space-time series estimation procedure}
\description{
\code{starma} fits a STARMA model to a space-time series. 
It is the central function for the estimation part of the three-stage iterative model building procedure.
}
\usage{
starma(data, wlist, ar, ma, iterate=1)

\method{print}{starma}(x, ...)
}
\arguments{
\item{data}{ a matrix or data frame containing the space-time series: row-wise should be the temporal observations, with each column corresponding to a site. }
\item{wlist}{ a list of the weight matrices for each k-th order neighbours, first one being the identity. }
\item{ar}{ either an integer specifying the maximum time lag of the AR part, or a matrix filled with 0 or 1 indicating whether 'row'-th tlag, 'col'-th slag AR parameter should be estimated. }
\item{ma}{ either an integer specifying the maximum time lag of the MA part, or a matrix filled with 0 or 1 indicating whether 'row'-th tlag, 'col'-th slag MA parameter should be estimated. }
\item{iterate}{ an integer specifying how many times the Kalman filter should be re-run on itself (see Details). }
\item{x}{ a \code{starma} class object. }
\item{...}{ unused }
}
\details{
The definition here used for STARMA models is the following:
\deqn{z_t = \sum_{k=1}^{p} \sum_{l=0}^{\lambda_k} \phi_{kl} W^{(l)} z_{t-k} + \sum_{k=1}^{q} \sum_{l=0}^{m_k} \theta_{kl} W^{(l)} \epsilon_{t-k} + \epsilon_t }{z(t) = phi[10]z(t-1) + phi[11]W[1]z(t-1) + ... + phi[pl]W[l]z(t-p) + theta[10]e(t-1) + theta[11]W[1]e(t-1) + ... + theta[qm]W[m]e(t-q) + e(t)}

\code{starma} uses a Kalman filter algorithm (Cipra and Motykova, 1987): the parameters are set as the state vector of the state space system, making the iterations of the algorithm estimate directly the parameters.
Thus, no optimization routine is required, making the algorithm extremely efficient time wise and computationally wise.
Furthermore, the code has been written in C++ using Rcpp and RcppArmadillo (Eddelbuettel and Sanderson, 2014).

Note that, as the residuals must be iteratively estimated when running the Kalman filter, a single run might lead to poor results when estimating an MA parameter.
Re-running the Kalman filter at least once, using the previously estimated parameters to add prior knowledge on the residuals leads to better estimates.
For STAR model (when no MA parameter needs be estimated), the function ignores the \code{iterate} argument.

One of the strength of this estimation function is that the user can to estimate as few parameters as needed, even at high time and/or space lags,
since the possibility to input a 1/0 matrix as AR and MA orders is given.
}
\value{
A list of class \code{starma} containing the following elements:
\item{phi}{ The estimated AR parameters }
\item{phi_sd}{ The corresponding standard errors }
\item{theta}{ The estimated MA parameters }
\item{theta_sd}{ The corresponding standard errors }
\item{sigma2}{ The white noise variance matrix estimated by the Kalman filter.
Note that, to achieve parcimony, only the mean of the diagonal elements should be kept (since the noise is supposed to be Gaussian anyway) }
\item{residuals}{ The estimated residuals of the model }
\item{loglik}{ The conditional log likelihood of the model }
\item{bic}{ The corresponding BIC }
\item{call}{ The function call }
\item{df}{ Degrees of freedom of the model: (nb of obs) - (nb of parameters) }
}
\author{Felix Cheysson}
\references{
- Cipra, T., & Motykova, I. (1987). Study on Kalman filter in time series analysis. Commentationes Mathematicae Universitatis Carolinae, 28(3).

- Dirk Eddelbuettel, Conrad Sanderson (2014). RcppArmadillo: Accelerating R with high-performance C++ linear algebra. Computational Statistics and Data Analysis, Volume 71, March 2014, pages 1054-1063. URL http://dx.doi.org/10.1016/j.csda.2013.02.005
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
stacf(sim, blist)
stpacf(sim, blist)

# Select parameters to estimate
ar <- matrix(0, 2, 2)
ar[ ,1] <- 1	# phi10 and phi20
ar[1,2] <- 1	# phi11
ma <- matrix(c(0,1), 1, 2)	# theta11

# Run the Kalman filter algorithm
model <- starma(sim, blist, ar, ma)
summary(model)
}
\keyword{starma}
\keyword{estimation}
