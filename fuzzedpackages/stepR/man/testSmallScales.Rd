\name{testSmallScales}
\alias{testSmallScales}
\alias{.testSmallScales}

\title{Test Small Scales}
\description{For developers only; users should look at the function \code{improveSmallScales} in the CRAN package \code{clampSeg}. Implements the second step of \acronym{HILDE} (\cite{Pein et al.}, 2020, Section III-B) in which an initial fit is tested for missed short events.}

\usage{
.testSmallScales(data, family, lengths = NULL, q, alpha, ...)
}

\arguments{
\item{data}{a numeric vector containing the observations}
\item{family}{a string specifying the assumed parametric family, currently \code{"LR"} and \code{"2Param"} are supported}
\item{lengths}{an integer vector giving the set of lengths, i.e. only intervals of these lengths will be considered. By default (\code{NULL}) \code{1:20} will be used for parametric family \code{"LR"} and \code{1:65} will be used for parametric family \code{"2Param"}}
\item{q}{either \code{NULL}, then the vector of critical values at level \code{alpha} will be computed from a Monte-Carlo simulation or a numeric vector giving the vector of critical values. Either \code{q} or \code{alpha} must be given. Otherwise, \code{alpha == 0.5} is chosen with a \code{\link{warning}}. This argument will be passed to \code{\link{critVal}} to obtain the needed critical values. Additional parameters for the computation of \code{q} can be specified in \code{\dots}, for more details see the documentation of \code{\link{critVal}}. Please note that by default the Monte-Carlo simulation will be saved in the workspace and on the file system, for more details see Section \cite{Storing of Monte-Carlo simulations} below}
\item{alpha}{a probability, i.e. a single numeric between 0 and 1, giving the significance level. Its choice is a trade-off between data fit and parsimony of the estimator. In other words, this argument balances the risks of missing change-points and detecting additional artefacts. For more details on this choice see (Frick et al., 2014, section 4) and (Pein et al., 2017, section 3.4). Either \code{q} or \code{alpha} must be given. Otherwise, \code{alpha == 0.5} is chosen with a \code{\link{warning}}}
\item{...}{there are two groups of further arguments:
\enumerate{
\item further parameters of the parametric family,
\item further parameters that will be passed to \code{\link{critVal}}. \code{\link{critVal}} will be called automatically with the number of observations \code{n = length(y)}, the arguments \code{family}, \code{intervalSystem}, \code{lengths}, \code{q} and \code{output} set. For these arguments no user interaction is required and possible, all other arguments of \code{\link{critVal}} can be passed additionally
}
}
}

\value{
a \code{\link{list}} with entries \code{jumps, addLeft, addRight, noDeconvolution, data, q}
}

\encoding{utf8}
\references{
Pein, F., Bartsch, A., Steinem, C., and Munk, A. (2020) Heterogeneous idealization of ion channel recordings - Open channel noise. Submitted.
}

\keyword{nonparametric}
