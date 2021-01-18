\name{BstrapTest}
\alias{BstrapTest}

\title{Bootstrap test for a single change-point}

\description{
Tests whether the underlying signal is smooth or contains at least one change-point. The smooth alternative is estimated by a (crossvalidated) kernel smoother. The single change-point alternative is estimated by \code{\link{estimateSingleCp}}. Its estimated jump size is used as a test statistic and the critical value is obtained by bootstrapping. More details can be found in the vignette.
}

\usage{
BstrapTest(y, bandwidth, nbandwidth = 30L, B = 500L, alpha = 0.05,
           kernel = c("epanechnikov", "gaussian", "rectangular",
                      "triangular", "biweight", "silverman"))
}

\arguments{
  \item{y}{a numeric vector containing the data points}
  \item{bandwidth}{the bandwidth, i.e. a numeric with values between \code{1 / length(y)} and \code{0.5}. If missing  \code{exp(seq(log(10 / length(y)), log(0.25), length.out = nbandwidth))} will be used. Crossvalidation will be performed if it is not a single numeric. Note that the test has almost no power when the bandwidth for the kernel smoother is too small, since then a change-point can be approximated well by a quickly changing smooth function.}
  \item{nbandwidth}{a single integer giving the number of bandwidths (see above) if \code{bandwidth} is missing}
   \item{B}{a single integer giving the number of bootstrap samples}
   \item{alpha}{a probability, i.e. a single numeric between 0 and 1, giving the significance level of the test}
  \item{kernel}{the kernel function, i.e. either a string or a function that takes a single numeric vector and returns the values of the kernel at those locations}
}

\value{
a \code{\link{list}} with the following components:\cr
- piecewiseSignal: the estimated signal with a single change-point\cr
- cp: the estimated change-point location\cr
- size: the estimated jump size\cr
- bandwidth: the selected bandwidth for the piecewise signal\cr
- bandwidthSmooth: the selected bandwidth for the smooth signal\cr
- smoothSignal: the estimated smooth signal\cr
- critVal: the by bootstrapping obtained critical value\cr
- pValue: the p-Value of the test\cr
- outcome: a boolean saying whether the test rejects the hypothesis of a smooth signal
}

\examples{
n <- 100
signal <- sin(2 * pi * 1:n / n)
signal[51:100] <- signal[51:100] + 5

y <- rnorm(n) + signal

# default bandwidth and kernel
test <- BstrapTest(y = y)

if (test$outcome) {
  # null hypothesis of a smooth signal is rejected
  estimatedSignal <- test$piecewiseSignal
} else {
  # null hypothesis of a smooth signal is accepted
  estimatedSignal <- test$smoothSignal
}

plot(y)
lines(signal)
lines(estimatedSignal, col = "red")

# fixed bandwidth
test <- BstrapTest(y = y, bandwidth = 0.1)

# user specified kernel
kernel <- function(x) 1 - abs(x) # triangular kernel
test <- BstrapTest(y = y, kernel = kernel)
}
  