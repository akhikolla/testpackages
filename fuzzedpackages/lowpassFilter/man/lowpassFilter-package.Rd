\name{lowpassFilter-package}
\alias{lowpassFilter-package}
\docType{package}

\title{Lowpass Filtering}

\description{Creates lowpass filters and offers further functionalities around them. Lowpass filters are commonly used in ion channel recordings.
}

\details{
The main function of this package is \code{\link{lowpassFilter}} which creates lowpass filters, currently only Bessel filters are supported. \code{\link{randomGeneration}} and \code{\link{randomGenerationMA}} allow to generate random numbers that are filtered, i.e. follow a model for ion channel recordings, see \cite{(Pein et al., 2018, 2020)}. \code{\link{getConvolution}}, \code{\link{getConvolutionJump}}, and \code{\link{getConvolutionPeak}} allow to compute the convolution of a signal with the kernel of a lowpass filter.
}

\encoding{utf8}
\references{
Pein, F., Bartsch, A., Steinem, C., and Munk, A. (2020) Heterogeneous idealization of ion channel recordings - Open channel noise. Submitted.

Pein, F., Tecuapetla-Gómez, I., \enc{Schütte}{Schuette}, O., Steinem, C., Munk, A. (2018) Fully-automatic multiresolution idealization for filtered ion channel recordings: flickering event detection. IEEE Trans. Nanobioscience, 17(3):300-320.

Pein, F. (2017) Heterogeneous Multiscale Change-Point Inference and its Application to Ion Channel Recordings.
PhD thesis, Georg-August-Universität Göttingen. http://hdl.handle.net/11858/00-1735-0000-002E-E34A-7.

Hotz, T., \enc{Schütte}{Schuette}, O., Sieling, H., Polupanow, T., Diederichsen, U., Steinem, C., and Munk, A. (2013) Idealizing
ion channel recordings by a jump segmentation multiresolution filter. IEEE Trans. Nanobioscience, 12(4):376-386.
}

\keyword{package}
\keyword{ts}
\keyword{nonparametric}

\seealso{\code{\link{lowpassFilter}}, \code{\link{randomGeneration}}, \code{\link{randomGenerationMA}}, \code{\link{getConvolution}}, \code{\link{getConvolutionJump}}, \code{\link{getConvolutionPeak}}}

\examples{
# creates a lowpass filter
filter <- lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1), sr = 1e4)
time <- 1:4000 / filter$sr

# creates a piecewise constant signal with a single peak
stepfun <- getSignalPeak(time, cp1 = 0.2, cp2 = 0.2 + 3 / filter$sr, 
                         value = 20, leftValue = 40, rightValue = 40)

# computes the convolution of the signal with the kernel of the lowpass filter
signal <- getConvolutionPeak(time, cp1 = 0.2, cp2 = 0.2 + 3 / filter$sr, 
                             value = 20, leftValue = 40, rightValue = 40,
                             filter = filter)

# generates random numbers that are filtered                             
data <- randomGenerationMA(n = 4000, filter = filter, signal = signal, noise = 1.4)

# generated data
plot(time, data, pch = 16)

# zoom into the single peak
plot(time, data, pch = 16, xlim = c(0.199, 0.202), ylim = c(19, 45))
lines(time, stepfun, col = "blue", type = "s", lwd = 2)
lines(time, signal, col = "red", lwd = 2)

# use of data randomGeneration instead
data <- randomGeneration(n = 4000, filter = filter, signal = signal, noise = 1.4)

# similar result
plot(time, data, pch = 16, xlim = c(0.199, 0.202), ylim = c(19, 45))
lines(time, stepfun, col = "blue", type = "s", lwd = 2)
lines(time, signal, col = "red", lwd = 2)
}
