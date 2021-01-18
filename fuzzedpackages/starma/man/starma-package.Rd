\name{starma-package}
\alias{starma-package}
\docType{package}
\title{
Space Time AutoRegressive Moving Average (STARMA) modelling for space-time series
}
\description{
This package aims to provide all the tools needed to identify, estimate and diagnose STARMA models for space-time series. 
It follows the three-stage iterative model building procedure developed by (Box and Jenkins, 1970) and extended to space-time modelling by (Pfeifer and Deutsch, 1980).
Designed with large datasets in mind, the package has been optimized by integrating C++ code via Rcpp and RcppArmadillo (Eddelbuettel and Sanderson, 2014).
Furthermore, the parameter estimation, which is usually computationally very expensive when using common optimization routines, uses a Kalman filter (see Cipra and Motykova, 1987), making it extremely efficient when dealing with large datasets.
}
\details{
\tabular{ll}{
Package: \tab starma\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2015-11-12\cr
License: \tab GPL-2\cr
}
The three stages of the iterative model building procedure are as follow, after centering the space-time series with \code{\link{stcenter}}:

- Identification: Using \code{\link{stacf}} and \code{\link{stpacf}}, the user should try to identify which parameters should be estimated. 

- Estimation: Use \code{\link{starma}} to estimate the parameters. 

- Diagnostic: Use \code{\link{stacf}}, \code{\link{stpacf}} and \code{\link{stcor.test}} to check whether the residuals of the models are similar to white noise. 

Refer to (Box and Jenkins, 1970) for details over the three-stage iterative model building procedure.
}
\author{
Felix Cheysson

Maintainer: Felix Cheysson <felix@cheysson.fr>
}
\references{
- Box, G. E. P., & Jenkins, G. M. (1970). Time Series Analysis: Forecasting and Control. Holden Day.

- Pfeifer, P., & Deutsch, S. (1980). A Three-Stage Iterative Procedure for Space-Time Modeling. Technometrics, 22(1), 35-47. doi:10.1080/00401706.1980.10486099

- Pfeifer, P., & Deutsch, S. (1981). Variance of the Sample Space-Time Autocorrelation Function. Journal of the Royal Statistical Society. Series B (Methodological), 43(1): 28-33.

- Cipra, T., & Motykova, I. (1987). Study on Kalman filter in time series analysis. Commentationes Mathematicae Universitatis Carolinae, 28(3).

- Dirk Eddelbuettel, Conrad Sanderson (2014). RcppArmadillo: Accelerating R with high-performance C++ linear algebra. Computational Statistics and Data Analysis, Volume 71, March 2014, pages 1054-1063. URL http://dx.doi.org/10.1016/j.csda.2013.02.005
}
\keyword{ package }
\keyword{ starma }
\keyword{ space-time }
\examples{
# Load spdep library to easily create weight matrices
library(spdep)

# Create a 5x5 regular grid which will be our lattice
sites <- matrix(0, 25, 2)
for (i in 1:5) {
	for (j in 1:5)
		sites[(i-1)*5 + j, ] <- c(i, j) - .5
}
plot(sites)

# Create a uniform first order neighbourhood
knb <- dnearneigh(sites, 0, 1)
plot(knb, sites)

# Lag the neighbourhood to create other order matrices
knb <- nblag(knb, 4)
klist <- list(order0=diag(25),
           order1=nb2mat(knb[[1]]),
           order2=nb2mat(knb[[2]]),
           order3=nb2mat(knb[[3]]),
           order4=nb2mat(knb[[4]]))
		   
# Simulate a STARMA(2;1) process
eps <- matrix(rnorm(200*25), 200, 25)
star <- eps
for (t in 3:200) {
	star[t,] <- (.4*klist[[1]] + .25*klist[[2]]) \%*\% star[t-1,] +
		(.25*klist[[1]]                ) \%*\% star[t-2,] +
		(            - .3*klist[[2]]) \%*\% eps[t-1,] +
		eps[t, ]
}

star <- star[101:200,]	# Remove first observations
star <- stcenter(star)	# Center and scale the dataset
										   
# Identify the process
stacf(star, klist)
stpacf(star, klist)

# Estimate the process
ar <- matrix(c(1, 1, 1, 0), 2, 2)
ma <- matrix(c(0, 1), 1, 2)
model <- starma(star, klist, ar, ma)
model
summary(model)

# Diagnose the process
stcor.test(model$residuals, klist, fitdf=4)
stacf(model$residuals, klist)
stpacf(model$residuals, klist)
}
