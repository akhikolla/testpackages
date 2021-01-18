\name{lowpassFilter}
\alias{lowpassFilter}
\alias{print.lowpassFilter}

\title{Lowpass filtering}

\description{
Creates a lowpass filter.
}

\usage{
lowpassFilter(type = c("bessel"), param, sr = 1, len = NULL, shift = 0.5)
\method{print}{lowpassFilter}(x, \dots)
}

\arguments{
\item{type}{a string specifying the type of the filter, currently only Bessel filters are supported}
\item{param}{a \code{\link{list}} specifying the parameters of the filter depending on \code{type}. For \code{"bessel"} the entries \code{pole} and \code{cutoff} have to be specified and no other named entries are allowed. \code{pole} has to be a single integer giving the number of poles (order). \code{cutoff} has to be a single positive numeric not larger than \code{1} giving the normalized cutoff frequency, i.e. the cutoff frequency (in the temporal domain) of the filter divided by the sampling rate}
\item{sr}{a single numeric giving the sampling rate}
\item{len}{a single integer giving the filter length of the truncated and digitised filter, see \cite{Value} for more details. By default (\code{NULL}) it is chosen such that the autocorrelation function is below \code{1e-3} at \code{len / sr} and at all lager lags \code{(len + i) / sr}, with \code{i} a positive integer}
\item{shift}{a single numeric between \code{0} and \code{1} giving a shift for the digitised filter, i.e. kernel and step are obtained by evaluating the corresponding functions at \code{(0:len + shift) / sr}}
\item{x}{the object}
\item{\dots}{for generic methods only}
}

\value{
An object of \code{\link{class}} \code{lowpassFilter}, i.e. a \code{\link{list}} that contains
\describe{
\item{\code{"type"}, \code{"param"}, \code{"sr"}, \code{"len"}}{the corresponding arguments}
\item{\code{"kernfun"}}{the kernel function of the filter, obtained as the Laplace transform of the corresponding transfer function}
\item{\code{"stepfun"}}{the step-response of the filter, i.e. the antiderivative of the filter kernel}
\item{\code{"acfun"}}{the autocorrelation function, i.e. the convolution of the filter kernel with itself}
\item{\code{"acAntiderivative"}}{the antiderivative of the autocorrelation function}
\item{\code{"truncatedKernfun"}}{the kernel function of the at \code{len / sr} truncated filter, i.e. \code{kernfun} truncated and rescaled such that the new kernel still integrates to \code{1}}
\item{\code{"truncatedStepfun"}}{the step-response of the at \code{len / sr} truncated filter, i.e. the antiderivative of the kernel of the truncated filter}
\item{\code{"truncatedAcfun"}}{the autocorrelation function of the at \code{len / sr} truncated filter, i.e. the convolution of the kernel of the truncated filter with itself}
\item{\code{"truncatedAcAntiderivative"}}{the antiderivative of the autocorrelation function of the at \code{len / sr} truncated filter}
\item{\code{"kern"}}{the digitised filter kernel normalised to one, i.e. \code{kernfun((0:len + shift) / sr) / sum(kernfun((0:len + shift) / sr))}}
\item{\code{"step"}}{the digitised step-response of the filter, i.e. \code{stepfun((0:len + shift) / sr)}}
\item{\code{"acf"}}{the discrete autocorrelation, i.e. \code{acfun(0:len / sr)}}
\item{\code{"jump"}}{the last index of the left half of the filter, i.e. \code{min(which(ret$step >= 0.5)) - 1L}, it indicates how much a jump is shifted in time by a convolution of the signal with the digitised kernel of the lowpassfilter; if all values are below \code{0.5}, \code{len} is returned with a \code{warning}}
\item{\code{"number"}}{for developers; an integer indicating the type of the filter}
\item{\code{"list"}}{for developers; a list containing precomputed quantities to recreate the filter in C++}
}
}

\seealso{\code{\link{filter}}}

\encoding{utf8}
\references{
Pein, F., Bartsch, A., Steinem, C., and Munk, A. (2020) Heterogeneous idealization of ion channel recordings - Open channel noise. Submitted.

Pein, F., Tecuapetla-Gómez, I., \enc{Schütte}{Schuette}, O., Steinem, C., Munk, A. (2018) Fully-automatic multiresolution idealization for filtered ion channel recordings: flickering event detection. IEEE Trans. Nanobioscience, 17(3):300-320.

Pein, F. (2017) Heterogeneous Multiscale Change-Point Inference and its Application to Ion Channel Recordings.
PhD thesis, Georg-August-Universität Göttingen. http://hdl.handle.net/11858/00-1735-0000-002E-E34A-7.

Hotz, T., \enc{Schütte}{Schuette}, O., Sieling, H., Polupanow, T., Diederichsen, U., Steinem, C., and Munk, A. (2013) Idealizing
ion channel recordings by a jump segmentation multiresolution filter. IEEE Trans. Nanobioscience, 12(4):376-386.
}

\examples{
filter <- lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 1e3 / 1e4),
                        sr = 1e4)

# filter kernel, truncated version
plot(filter$kernfun, xlim = c(0, 20 / filter$sr))
t <- seq(0, 20 / filter$sr, 0.01 / filter$sr)
# truncated version looks very similar
lines(t, filter$truncatedKernfun(t), col = "red")

# filter$len (== 11) is chosen automatically
# this ensures that filter$acf < 1e-3 for this lag and at all larger lags
plot(filter$acfun, xlim = c(0, 20 / filter$sr), ylim = c(-0.003, 0.003))
abline(h = 0.001, lty = "22")
abline(h = -0.001, lty = "22")

abline(v = (filter$len - 1L) / filter$sr, col = "grey")
abline(v = filter$len / filter$sr, col = "red")

# filter with sr == 1
filter <- lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 1e3 / 1e4))

# filter kernel and its truncated version
plot(filter$kernfun, xlim = c(0, 20 / filter$sr))
t <- seq(0, 20 / filter$sr, 0.01 / filter$sr)
# truncated version looks very similar
lines(t, filter$truncatedKernfun(t), col = "red")
# digitised filter
points((0:filter$len + 0.5) / filter$sr, filter$kern, col = "red", pch = 16)

# without a shift
filter <- lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 1e3 / 1e4),
                        shift = 0)
# filter$kern starts with zero
points(0:filter$len / filter$sr, filter$kern, col = "blue", pch = 16)

# much shorter filter
filter <- lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 1e3 / 1e4),
                        len = 4L)
points((0:filter$len + 0.5) / filter$sr, filter$kern, col = "darkgreen", pch = 16)
}

\keyword{ts}
