\name{convolve}
\alias{convolve}
\alias{.convolve}

\title{Time discrete convolution}
\description{For developers only; computes a time discrete convolution.}

\usage{
.convolve(val, kern)
}

\arguments{
\item{val}{a numeric vector giving the values}
\item{kern}{a numeric vector giving the time discrete kernel}
}

\value{
A numeric vector giving the convolution.
}

\seealso{\code{\link{lowpassFilter}}}
