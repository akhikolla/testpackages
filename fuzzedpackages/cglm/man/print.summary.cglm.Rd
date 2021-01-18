\name{print.summary.cglm}

\alias{print.summary.cglm}

\title{
Prints summary of instrumental variable estimation
}

\description{
This is a \code{print} method for class \code{"summary.cglm"}.
}

\usage{
\method{print}{summary.cglm}(x, digits=max(3L, getOption("digits")-3L), 
  signif.stars=getOption("show.signif.stars"), ...)
}

\arguments{
  \item{x}{
an object of class \code{"summary.cglm"}.
}
  \item{digits}{
the number of significant digits to use.
}
  \item{signif.stars}{
logical. If TRUE, "significance stars" are printed for each coefficient.
}
  \item{\dots}{
not used.
}
}

\author{
Arvid Sjolander
}


\examples{

##See documentation for cglm.

}
