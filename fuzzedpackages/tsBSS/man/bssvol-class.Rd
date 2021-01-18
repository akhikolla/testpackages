\name{bssvol}
\docType{class}
\alias{bssvol}
\alias{print.bssvol}
\alias{plot.bssvol}
\title{
Class: bssvol
}

\description{
Class bssvol (blind source separation in stochastic volatility processes) with methods print.bssvol (prints an object of class bssvol) and plot.bss (plots an object of class bssvol).

Class also inherits methods from the class bss in package \code{\link[JADE]{JADE}}: for extracting the components of an object of class bssvol (\code{\link[JADE:bss.components]{bss.components}}) and the coefficients of an object of class bssvol (\code{\link[JADE:coef.bss]{coef.bss}}). 
}

\usage{
\method{print}{bssvol}(x, ...)

\method{plot}{bssvol}(x, ...)
}

\arguments{
  \item{x}{An object of class bssvol.}
  \item{\dots}{Further arguments to be passed to or from methods.}
}

\author{
Markus Matilainen
}

\seealso{
\code{\link[JADE]{JADE}}, \code{\link[JADE:bss.components]{bss.components}}, \code{\link[JADE:coef.bss]{coef.bss}}
}

\keyword{methods}