\name{spray-class}
\docType{class}
\alias{spray-class}
\title{Class \dQuote{spray}}
\description{The formal S4 class for sprays.}
\section{Objects from the Class}{

Objects \emph{can} be created by calls of the form \code{new("spray",
  ...)} but this is not encouraged.  Use functions \code{spray()} or
 \code{as.spray()} instead.
}
\section{Slots}{
  \describe{
    \item{\code{index}:}{Index matrix}
    \item{\code{value}:}{Numeric vector holding coefficients}
  }
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{spray}}}
\keyword{classes}
