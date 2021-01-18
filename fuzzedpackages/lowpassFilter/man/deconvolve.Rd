\name{deconvolve}
\alias{deconvolve}
\alias{.deconvolveJump}
\alias{.deconvolvePeak}
\alias{deconvolveJump}
\alias{deconvolvePeak}

\title{Deconvolution of a single jump / isolated peak}
\description{For developers only; computes the deconvolution of a single jump or an isolated peak assuming that the observations are lowpass filtered. More details are given in \cite{(Pein et al., 2018)}.}

\usage{
.deconvolveJump(grid, observations, time, leftValue, rightValue,
                typeFilter, inputFilter, covariances) 
.deconvolvePeak(gridLeft, gridRight, observations, time, leftValue, rightValue,
                typeFilter, inputFilter, covariances, tolerance)
}

\arguments{
\item{grid, gridLeft, gridRight}{numeric vectors giving the potential time points of the single jump, of the left and right jump points of the peak, respectively}
\item{observations}{a numeric vector giving the observed data}
\item{time}{a numeric vector of length \code{length(observations)} giving the time points at which the \code{observations} are observed}
\item{leftValue, rightValue}{single numerics giving the value (conductance level) before and after the jump / peak, respectively}
\item{typeFilter, inputFilter}{a description of the assumed lowpass filter, usually computed by \code{\link{lowpassFilter}}}
\item{covariances}{a numeric vector giving the (regularized) covariances of the observations}
\item{tolerance}{a single numeric giving a tolerance for the decision whether the left jump point is smaller than the right jump point}
}

\value{
For \code{.deconvolveJump} a single numeric giving the jump point. For \code{.deconvolvePeak} a list containing the entries \code{left}, \code{right} and \code{value} giving the left and right jump point and the value of the peak, respectively.
}

\seealso{\code{\link{lowpassFilter}}}

\encoding{utf8}
\references{
Pein, F., Tecuapetla-Gómez, I., \enc{Schütte}{Schuette}, O., Steinem, C., Munk, A. (2018) Fully-automatic multiresolution idealization for filtered ion channel recordings: flickering event detection. IEEE Trans. Nanobioscience, 17(3):300-320.
}
