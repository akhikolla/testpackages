\name{helpFunctionsFilter}
\alias{getConvolution}
\alias{getSignalJump}
\alias{getConvolutionJump}
\alias{getSignalPeak}
\alias{getConvolutionPeak}


\title{Convolved piecewise constant signals}
\description{Creates piecewise constant signals with a single jump / peak. Computes the convolution of piecewise constant signals with the kernel of a lowpass filter.}

\usage{
getConvolution(t, stepfun, filter, truncated = TRUE)
getSignalJump(t, cp, leftValue, rightValue)
getConvolutionJump(t, cp, leftValue, rightValue, filter, truncated = TRUE)
getSignalPeak(t, cp1, cp2, value, leftValue, rightValue) 
getConvolutionPeak(t, cp1, cp2, value, leftValue, rightValue, filter, truncated = TRUE) 
}

\arguments{
\item{t}{a numeric vector giving the time points at which the signal / convolution should be computed}
\item{stepfun}{specification of the piecewise constant signal, i.e. a \code{\link{data.frame}} with named arguments \code{leftEnd}, \code{rightEnd} and \code{value} giving the start and end points of the constant segments and the values on the segments, for instance an object of class \code{stepblock} as available by the package 'stepR'}
\item{cp, cp1, cp2}{a single numeric giving the location of the single, first and second jump point, respectively}
\item{value, leftValue, rightValue}{a single numeric giving the function value at, before and after the peak / jump, respectively}
\item{filter}{an object of class \code{\link{lowpassFilter}} giving the analogue lowpass filter}
\item{truncated}{a single logical (not NA) indicating whether the signal should be convolved with the truncated or the untruncated filter kernel}
}

\value{
a numeric of length \code{length(t)} giving the signal / convolution at time points \code{t}
}

\seealso{\code{\link{lowpassFilter}}}

\encoding{utf8}
\references{
Pein, F., Bartsch, A., Steinem, C., and Munk, A. (2020) Heterogeneous idealization of ion channel recordings - Open channel noise. Submitted.

Pein, F., Tecuapetla-Gómez, I., \enc{Schütte}{Schuette}, O., Steinem, C., Munk, A. (2018) Fully-automatic multiresolution idealization for filtered ion channel recordings: flickering event detection. IEEE Trans. Nanobioscience, 17(3):300-320.

Pein, F. (2017) Heterogeneous Multiscale Change-Point Inference and its Application to Ion Channel Recordings.
PhD thesis, Georg-August-Universität Göttingen. http://hdl.handle.net/11858/00-1735-0000-002E-E34A-7.
}

\examples{
# creating and plotting a signal with a single jump at 0 from 0 to 1
time <- seq(-2, 13, 0.01)
signal <- getSignalJump(time, 0, 0, 1)
plot(time, signal, type = "l")

# setting up the filter
filter <- lowpassFilter(param = list(pole = 4, cutoff = 0.1))

# convolution with the truncated filter
convolution <- getConvolutionJump(time, 0, 0, 1, filter)
lines(time, convolution, col = "red")

# without truncating the filter, looks almost equal
convolution <- getConvolutionJump(time, 0, 0, 1, filter, truncated = FALSE)
lines(time, convolution, col = "blue")


# creating and plotting a signal with a single peak with jumps
# at 0 and at 3 from 0 to 1 to 0
time <- seq(-2, 16, 0.01)
signal <- getSignalPeak(time, 0, 3, 1, 0, 0)
plot(time, signal, type = "l")

# convolution with the truncated filter
convolution <- getConvolutionPeak(time, 0, 3, 1, 0, 0, filter)
lines(time, convolution, col = "red")

# without truncating the filter, looks almost equal
convolution <- getConvolutionPeak(time, 0, 3, 1, 0, 0, filter, truncated = FALSE)
lines(time, convolution, col = "blue")


# doing the same with getConvolution
# signal can also be an object of class stepblock instead,
# e.g. constructed by stepR::stepblock
signal <- data.frame(value = c(0, 1, 0), leftEnd = c(-2, 0, 3), rightEnd = c(0, 3, 16))

convolution <- getConvolution(time, signal, filter)
lines(time, convolution, col = "red")

convolution <- getConvolution(time, signal, filter, truncated = FALSE)
lines(time, convolution, col = "blue")


# more complicated signal
time <- seq(-2, 21, 0.01)
signal <- data.frame(value = c(0, 10, 0, 50, 0), leftEnd = c(-2, 0, 3, 6, 8),
                     rightEnd = c(0, 3, 6, 8, 21))

convolution <- getConvolution(time, signal, filter)
plot(time, convolution, col = "red", type = "l")

convolution <- getConvolution(time, signal, filter, truncated = FALSE)
lines(time, convolution, col = "blue")
}

\keyword{nonparametric}
