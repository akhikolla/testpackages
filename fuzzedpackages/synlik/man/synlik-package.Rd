\name{synlik-package}
\alias{synlik-package}
\docType{package}
\title{
Synthetic Likelihood Methods for Intractable Likelihoods
}
\description{
Package that provides Synthetic Likelihood methods for intractable likelihoods. The package is meant to be as general purpose
as possible: as long as you are able to simulate data from your model you should be able to fit it.
}
\details{
\tabular{ll}{
Package: \tab synlik\cr
Type: \tab Package\cr
Version: \tab 0.1.2 \cr
Date: \tab 2018-05-22\cr
License: \tab GPL (>=2)\cr
}
The package allows users to create objects of class \code{synlik} (S4), which are essentially constituted of a \code{simulator} function and 
a function (\code{summaries}) that transforms the data into summary statistics. The \code{simulator} can output any kind of data (vector, list, etc)
and this will be passed directly to the \code{summaries} function. This allow the package to fit a large variety of models. \cr

Once the model of interest has been set up as a \code{synlik} object, it is possible several methods on it. The function most useful function is \code{slik}, which can be used to evaluate the synthetic likelihood. The \code{slice.synlik} function allows to obtain and plot slices of the synthetic likelihood with respect to model parameters. It is possible to simulate data or statistics from the model using the generic \code{simulate}, and to check the normality of the statistics using the \code{checkNorm} function. Unknow parameters can be estimated by MCMC, through the \code{smcmc} function. This function will return an object of class \code{smcmc} (S4), which contains all the inputs and results of the MCMC procedure. \cr

Many functions in the package support parallel simulation on multiple cores.
}
\author{
Matteo Fasiolo and Simon N. Wood

Maintainer: Matteo Fasiolo <matteo.fasiolo@gmail.com>
}
\references{
Simon N Wood. Statistical inference for noisy nonlinear ecological dynamic systems. Nature, 466(7310):1102--1104, 2010.
}
\keyword{ Intractable Likelihood, Simulation-based inference. }
\seealso{
For some examples see the Vignettes (type \code{vignette("synlik")}).
}
\examples{

\dontrun{
#### Here I put a simple example, 
#### if you want to see more type: vignette("synlik")
}

#### Create synlik object
ricker_sl <- synlik(simulator = rickerSimul,
                    summaries = rickerStats,
                    param = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)),
                    extraArgs = list("nObs" = 50, "nBurn" = 50),
                    plotFun = function(input, ...){
                               plot(drop(input), type = 'l', ylab = "Pop", xlab = "Time", ...)
                            }
)

#### Simulate from the object
ricker_sl@data <- simulate(ricker_sl)
ricker_sl@extraArgs$obsData <- ricker_sl@data

#### Simulate statistics (each row is a vector of statistics)
simulate(ricker_sl, seed = 523, nsim = 10, stats = TRUE)

#### Plotting the data
plot(ricker_sl)

#### Checking multivariate normality of the statistics
checkNorm(ricker_sl)

#### Evaluate the likelihood
set.seed(4234)
slik(ricker_sl, 
     param  = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)),
     nsim    = 1e3)

#### Plotting a slice of the log-Likelihood possibly using multiple cores
slice(object = ricker_sl, 
      ranges = list("logR" = seq(3.5, 3.9, by = 0.02),
                    "logPhi" = seq(2, 2.6, by = 0.02),
                    "logSigma" = seq(-2, -0.5, by = 0.05)), 
      param = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)), 
      nsim = 500, multicore = FALSE)

#### MCMC estimation possibly using multiple cores
set.seed(4235)
ricker_sl <- smcmc(ricker_sl, 
                   initPar = c(3.2, -1, 2.6),
                   niter = 50, 
                   burn = 3,
                   priorFun = function(input, ...) 0, 
                   propCov = diag(c(0.1, 0.1, 0.1))^2, 
                   nsim = 1e3, 
                   multicore = FALSE)

# Continue with additional 50 iterations
ricker_sl <- continue(ricker_sl, niter = 50)

# Plotting results on transformed scale (exponential)
trans <- rep("exp", 3)
names(trans) <- names(ricker_sl@param)

plot(ricker_sl)
}
