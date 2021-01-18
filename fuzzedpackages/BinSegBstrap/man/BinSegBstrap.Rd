\name{BinSegBstrap}
\alias{BinSegBstrap}

\title{Estimates a piecewise smooth signal}

\description{
A piecewise smooth signal is estimated by applying \code{\link{BstrapTest}} recursively (binary segmentation approach). The final estimator is estimated by kernel smoothing on each segment separately; a joint bandwidth is selected by crossvalidation. More details can be found in the vignette.
}

\usage{
BinSegBstrap(y, bandwidth, nbandwidth = 30L, B = 500L, alpha = 0.05,
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
- est: the estimated signal\cr
- cps: the estimated change-point locations\cr
- bandwidth: the selected bandwidth
}

\examples{
n <- 200
signal <- sin(2 * pi * 1:n / n)
signal[51:100] <- signal[51:100] + 5
signal[151:200] <- signal[151:200] + 5

y <- rnorm(n) + signal

# default bandwidth and kernel
est <- BinSegBstrap(y = y)

plot(y)
lines(signal)
lines(est$est, col = "red")

# fixed bandwidth
est <- BinSegBstrap(y = y, bandwidth = 0.1)

# user specified kernel
kernel <- function(x) 1 - abs(x) # triangular kernel
est <- BinSegBstrap(y = y, kernel = kernel)
}
  