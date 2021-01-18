\name{scrBesag}
\alias{scrBesag}
\title{Simultaneous credible band computation (Besag algorithm)}
\usage{scrBesag(samples, level=0.95)}
\description{Simultaneous credible band computation}
\details{Calculates a series of simultaneous credible bounds for one parameter type.}
\value{matrix with \sQuote{lower} and \sQuote{upper} rows}
\references{J. Besag, P. Green, D. Higdon, K. Mengersen (1995):
Bayesian Methods for Nonlinear Classification and Regression,
\emph{Statistical Science} \bold{10}/1, 3--66,
\doi{10.1214/ss/1177010123}}
\author{Thomas Kneib}
\keyword{htest}
\arguments{\item{samples}{m by n matrix where m is the number of parameters,
n is the number of samples and  hence each (multivariate) sample is a column
in the matrix \code{samples}}
\item{level}{the credible level (default: 0.95)}}
