\name{derivatives.est}
\alias{derivatives.est}
\title{Estimation of the First Two Derivatives for Functional Data}
\usage{
derivatives.est(dataf, range = NULL, d = 101, spar = NULL, deriv = c(0,1))
}
\arguments{
\item{dataf}{Functional dataset, represented by a \code{dataf} object of their arguments
and functional values. \code{m} stands for the number of functions.}

\item{range}{The common range of the domain where the functions \code{dataf} are observed.
Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
\code{dataf}.}

\item{d}{Grid size to which all the functional data are transformed. For computation, 
all functional observations are first transformed into vectors of their functional values of length \code{d}
corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
points are reconstructed using linear interpolation, and extrapolation.}

\item{spar}{If provided, this parameter is passed to functions \code{D1ss} and \code{D2ss} from package \code{sfsmisc}
as the value of the smoothing spline parameter in order to numerically approximate
the derivatives of \code{dataf}.}

\item{deriv}{A vector composed of \code{0}, \code{1}, and \code{2} of the demanded 
functional values / derivatives of the functions in the rows of \code{dataf}.
\code{0} stands for the functional values, \code{1} for the first derivatives, 
\code{2} for the second derivatives.}
}
\value{
A multivariate \code{dataf} object of the functional values and / or the derivatives of \code{dataf}. 
The dimensionality of the vector-valued functional data is \code{nd}. The arguments of the data are all equal to 
an equi-distant grid of \code{d} points in the domain given by \code{range}. \code{nd} is the demanded number 
of derivatives at the output, i.e. the length of the vector \code{deriv}.
}
\description{
Returns the estimated values of derivatives of functional data.
}
\details{
If the input \code{dataf} is a functional random sample of size \code{m}, 
the function returns a \code{dataf} object of \code{nd}-dimensional functional data, where 
in the elements of the vector-valued functional data represent the estimated values of the 
derivatives of \code{dataf}. All derivatives are evaluated at an equi-distant grid of \code{d}
points in the domain given by \code{range}. \code{nd} here stands for \code{1}, \code{2} or \code{3}, 
depending on how many derivatives of \code{dataf} are
requested to be computed. For the estimation, functions \code{D1ss} and \code{D2ss} from the package
\code{sfsmisc} are utilized.
}
\examples{
dataf = dataf.population()$dataf
derivatives.est(dataf,deriv=c(0,1,2))
}
\seealso{
\code{\link[sfsmisc:D2ss]{D1ss}} in package \pkg{sfsmisc}

\code{\link[sfsmisc:D2ss]{D2ss}} in package \pkg{sfsmisc}
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{derivatives}
\keyword{functional}
\keyword{kernel}
