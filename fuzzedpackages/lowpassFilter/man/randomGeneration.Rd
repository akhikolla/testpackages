\name{randomGeneration}
\alias{randomGeneration}
\alias{randomGenerationMA}


\title{Random number generation}
\description{Generate random numbers that are filtered. Both, signal and noise, are convolved with the given lowpass filter, see details. Can be used to generate synthetic data resembling ion channel recordings, please see \cite{(Pein et al., 2018, 2020)} for the exact models.}

\usage{
randomGeneration(n, filter, signal = 0, noise = 1, oversampling = 100L, seed = n,
                 startTime = 0, truncated = TRUE)
randomGenerationMA(n, filter, signal = 0, noise = 1, seed = n,
                   startTime = 0, truncated = TRUE)                 
}

\arguments{
\item{n}{a single positive integer giving the number of observations that should be generated}
\item{filter}{an object of class \code{\link{lowpassFilter}} giving the analogue lowpass filter}
\item{signal}{either a numeric of length 1 or of length \code{n} giving the convolved signal, i.e. the mean of the random numbers, or an object that can be passed to \code{getConvolution}, i.e. an object of class \code{stepblock}, see \cite{Examples}, giving the signal that will be convolved with the kernel of the lowpass filter \code{filter}}
\item{noise}{for \code{randomGenerationMA} a single positive finite numeric giving the constant noise level, for \code{randomGeneration} either a numeric of length 1 or of length \code{(n + filter$len - 1L) * oversampling} or an object of class \code{stepblock}, see \cite{Examples}, giving the noise of the random errors, see \cite{Details}}
\item{oversampling}{a single positive integer giving the factor by which the errors should be oversampled, see \cite{Details}}
\item{seed}{will be passed to \code{\link{set.seed}} to set a seed, \code{\link{set.seed}} will not be called if this argument is set to \code{"no"}, i.e. a single value, interpreted as an \code{\link{integer}}, \code{NULL} or \code{"no"}}
\item{startTime}{a single finite numeric giving the time at which sampling should start}
\item{truncated}{a single logical (not NA) indicating whether the signal should be convolved with the truncated or the untruncated filter kernel}
}

\value{
a numeric vector of length \code{n} giving the generated random numbers
}

\details{
As discussed in \cite{(Pein et al., 2018)} and \cite{(Pein et al., 2020)}, in ion channel recordings the recorded data points can be modelled as equidistant sampled at rate \code{filter$sr} from the convolution of a piecewise constant signal perturbed by Gaussian white noise scaled by the noise level with the kernel of an analogue lowpass filter. The noise level is either constant (homogeneous noise, see \cite{(Pein et al., 2018)}) or itself varying (heterogeneous noise, see \cite{(Pein et al., 2020)}). \code{randomGeneration} and \code{randomGenerationMA} generate synthetic data from such models. \code{randomGeneration} allows homogeneous and heterogeneous noise, while \code{randomGenerationMA} only allows homogeneous noise, i.e. \code{noise} has to be a single numeric giving the constant noise level. The resulting observations represent the conductance at time points \code{startTime + 1:n / filter$sr}.

The generated observations are the sum of a convolved signal evaluated at those time points plus centred Gaussian errors that are correlated (coloured noise), because of the filtering, and scaled by the noise level. The convolved signal evaluated at those time points can either by specified in \code{signal} directly or \code{signal} can specify a piecewise constant signal that will be convolved with the \code{filter} using \code{\link{getConvolution}} and evaluated at those time points. \code{randomGenerationMA} computes a moving average process with the desired autocorrelation to generate random errors. \code{randomGeneration} oversamples the error, i.e. generates errors at time points \code{startTime + (seq(1 - filter$len + 1 / oversampling, n, 1 / oversampling) - 1 / 2 / oversampling) / filter$sr}, which will then be convolved with the \code{filter}. For this function \code{noise} can either give the noise levels at those oversampled time points or specify a piecewise constant function that will be automatically evaluated at those time points.
}

\seealso{\code{\link{lowpassFilter}}, \code{\link{getConvolution}}}

\encoding{utf8}
\references{
Pein, F., Bartsch, A., Steinem, C., and Munk, A. (2020) Heterogeneous idealization of ion channel recordings - Open channel noise. Submitted.

Pein, F., Tecuapetla-Gómez, I., \enc{Schütte}{Schuette}, O., Steinem, C., Munk, A. (2018) Fully-automatic multiresolution idealization for filtered ion channel recordings: flickering event detection. IEEE Trans. Nanobioscience, 17(3):300-320.

Pein, F. (2017) Heterogeneous Multiscale Change-Point Inference and its Application to Ion Channel Recordings.
PhD thesis, Georg-August-Universität Göttingen. http://hdl.handle.net/11858/00-1735-0000-002E-E34A-7.
}

\examples{
filter <- lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1), sr = 1e4)
time <- 1:4000 / filter$sr
stepfun <- getSignalPeak(time, cp1 = 0.2, cp2 = 0.2 + 3 / filter$sr, 
                         value = 20, leftValue = 40, rightValue = 40)
signal <- getConvolutionPeak(time, cp1 = 0.2, cp2 = 0.2 + 3 / filter$sr, 
                             value = 20, leftValue = 40, rightValue = 40, filter = filter)
data <- randomGenerationMA(n = 4000, filter = filter, signal = signal, noise = 1.4)

# generated data
plot(time, data, pch = 16)

# zoom into the single peak
plot(time, data, pch = 16, xlim = c(0.199, 0.202), ylim = c(19, 45))
lines(time, stepfun, col = "blue", type = "s", lwd = 2)
lines(time, signal, col = "red", lwd = 2)

# use of randomGeneration instead
data <- randomGeneration(n = 4000, filter = filter, signal = signal, noise = 1.4)

# similar result
plot(time, data, pch = 16, xlim = c(0.199, 0.202), ylim = c(19, 45))
lines(time, stepfun, col = "blue", type = "s", lwd = 2)
lines(time, signal, col = "red", lwd = 2)

## heterogeneous noise
# manual creation of an object of class 'stepblock'
# instead the function stepblock in the package stepR can be used
noise <- data.frame(leftEnd = c(0, 0.2, 0.2 + 3 / filter$sr),
                    rightEnd = c(0.2, 0.2 + 3 / filter$sr, 0.4),
                    value = c(1, 30, 1))
attr(noise, "x0") <- 0
class(noise) <- c("stepblock", class(noise))

data <- randomGeneration(n = 4000, filter = filter, signal = signal, noise = noise)

plot(time, data, pch = 16, xlim = c(0.199, 0.202), ylim = c(19, 45))
lines(time, stepfun, col = "blue", type = "s", lwd = 2)
lines(time, signal, col = "red", lwd = 2)
}

\keyword{nonparametric}
