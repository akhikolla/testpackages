\name{DiffusionRjgqd-package}
\alias{DiffusionRjgqd-package}
\alias{DiffusionRjgqd}
\docType{package}
\title{
A package for Performing Inference and Analysis on Generalized Quadratic Jump Diffusion Processes (JGQDs).
}
\description{
\bold{DiffusionRjgqd} is a toolbox for performing analysis and inference on a class of jump diffusion processes with quadratic drift and diffusion with state-dependent and time inhomogeneous jump mechanisms. The package consists of functions for performing likelihood based inference and transitional density approximations for both 1D and 2D JGQDs.
}
\details{
\tabular{ll}{
Package: \tab DiffusionRjgqd\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2015-12-01\cr
License: \tab GPL (>= 2)\cr
}
The package is designed around an interface whereby the user supplies standard R-functions dictating the functional form of the coefficients of the JGQD after which the package handles all the necessary mathematics and algorithmic construction. Furthermore, computational efficiency is optimized by constructing algorithms in C++ using the \bold{Rcpp} and \bold{RcppArmadillo} libraries.

Functions included in the package:
\tabular{lcl}{
\code{\link{BiJGQD.density}} \tab:\tab Generate the transitional density of a 2D JGQD.\cr
\code{\link{BiJGQD.mcmc}}*   \tab:\tab Conduct inference via MCMC on a 2D JGQD.\cr
\code{\link{JGQD.density}}   \tab:\tab Generate the transitional density of a 1D JGQD.\cr
\code{\link{JGQD.mcmc}}*     \tab:\tab Conduct inference via MCMC on a 1D JGQD.\cr
}
* Functions use C++.
}

\author{
Etienne A.D. Pienaar: \email{etiennead@gmail.com}

Maintainer: Etienne A.D. Pienaar (\email{etiennead@gmail.com})
}
\references{
Updates available on GitHub at \url{https://github.com/eta21}.
}
\keyword{ package }
\keyword{ C++}
\seealso{
\code{\link{BiJGQD.mcmc}}, \code{\link{JGQD.mcmc}}, \code{\link{JGQD.dic}},\code{\link{JGQD.remove}} and \code{\link{JGQD.density}}.
}

\examples{
\dontrun{
example(JGQD.density)
example(BiJGQD.density)
}
}
