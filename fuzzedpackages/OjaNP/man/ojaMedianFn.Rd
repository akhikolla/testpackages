\name{ojaMedianFn}
\alias{ojaMedianFn}
\title{ Value of the Oja Criterion}
\description{
 The function returns the value of the Oja criterion function for a given point.
}
\usage{
ojaMedianFn(X, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data set or data matrix.}
  \item{x}{a numeric vector with the coordinates of interest.}
}

\value{
a numeric value containing the ciriterion evaluated at \code{x} for the data \code{X}.
}
\references{
\cite{Oja, H. (1983), Descriptive statistics for multivariate distributions, \emph{Statistics and Probability Letters}, \bold{1}, 327--332. }

\cite{Ronkainen, T., Oja, H. and Orponen, P. (2002), Computation of the multivariate Oja median, in Dutter R.,  Filzmoser P.,Gather U. and Rousseeuw, P. J.: \emph{Developments in Robust Statistics}, Heidelberg: Springer, 344--359. }
}

\author{ported to R by Daniel Fischer. Original C++ code by Tommi Ronkainen.}

\seealso{\code{\link{ojaMedian}} }
\examples{
data(biochem)
X <- as.matrix(biochem[,1:2])
x <- ojaMedian(X)
ojaMedianFn(X, x)
ojaMedianFn(X, c(1.1, 0.4))
}
\keyword{multivariate}
\keyword{nonparametric}
