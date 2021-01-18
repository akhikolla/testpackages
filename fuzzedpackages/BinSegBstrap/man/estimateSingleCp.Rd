\name{estimateSingleCp}
\alias{estimateSingleCp}

\title{Estimation of a single change-point}

\description{
Estimates a single change-point in an otherwise smooth function. The change-point location is estimated as the maximum of the differences of left and right sided running means. The estimate left and right of the change-point are obtained by kernel smoothers. Windows of the running mean and kernel bandwidth are chosen by crossvalidation. More details can be found in the vignette.
}

\usage{
estimateSingleCp(y, bandwidth, nbandwidth = 30L,
                 kernel = c("epanechnikov", "gaussian", "rectangular",
                            "triangular", "biweight", "silverman"))
}

\arguments{
  \item{y}{a numeric vector containing the data points}
  \item{bandwidth}{the bandwidth, i.e. a numeric with values between \code{1 / length(y)} and \code{0.5}. If missing  \code{exp(seq(log(2 / length(y)), log(0.25), length.out = nbandwidth))} will be used. Crossvalidation will be performed if it is not a single numeric}
  \item{nbandwidth}{a single integer giving the number of bandwidths (see above) if \code{bandwidth} is missing}
  \item{kernel}{the kernel function, i.e. either a string or a function that takes a single numeric vector and returns the values of the kernel at those locations}
}

\value{
a \code{\link{list}} with the following components:\cr
- est: the estimated function with a single change-point\cr
- cp: the estimated change-point location\cr
- size: the estimated jump size\cr
- bandwidth: the selected bandwidth
}

\examples{
n <- 100
signal <- sin(2 * pi * 1:n / n)
signal[51:100] <- signal[51:100] + 5

y <- rnorm(n) + signal

# default bandwidth and kernel
est <- estimateSingleCp(y = y)

plot(y)
lines(signal)
lines(est$est, col = "red")

# fixed bandwidth
est <- estimateSingleCp(y = y, bandwidth = 0.1)

# user specified kernel
kernel <- function(x) 1 - abs(x) # triangular kernel
est <- estimateSingleCp(y = y, kernel = kernel)
}
  