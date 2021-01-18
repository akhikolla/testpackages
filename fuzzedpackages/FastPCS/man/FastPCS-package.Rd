\name{FastPCS-package}
\docType{package}
\title{FastPCS estimator of location and scatter.}
\description{
Uses the FastPCS algorithm to compute the robust PCS estimator of location and scatter.
}
\details{
\tabular{ll}{
Package: \tab FastPCS\cr
Type: \tab Package\cr
Version: \tab 0.0.2\cr
Date: \tab 2013-01-13\cr
Suggests: \tab mvtnorm\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

Index:
\preformatted{
FastPCS                 Function to compute the robust FastPCS estimator of location and scatter.
FPCSnumStarts               Internal function used to compute the number of starting points used by FastPCS.
quanf                   Internal function used to compute h, the minimum number of observations expected to be uncontaminated.
plot.FastPCS		Plots the robust distances outputted by a FastPCS model.
}
}
\references{
	Vakili, K. and Schmitt, E. (2014).
	Finding multivariate outliers with FastPCS.
	\emph{Computational Statistics \& Data Analysis}.
  	Vol. 69, pp 54--66.
	(http://arxiv.org/abs/1301.2053)
}
\author{
Kaveh Vakili [aut, cre],
Maintainer: Kaveh Vakili <vakili.kaveh.email@gmail.com>
}
\keyword{package}
