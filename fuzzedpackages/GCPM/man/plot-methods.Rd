\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot}
\alias{plot,ANY-method}
\alias{plot,GCPM-method}
\title{Plot of the Portfolio Loss Distribution}
\description{
Plot of the estimated pdf of the portfolio loss distribution.
}
\usage{
plot(x,y,...)
}
\arguments{
  \item{x}{Object of class \code{GCPM}}
  \item{y}{plot unit for losses (x-axis), default value = 1}
  \item{...}{Further arguments such as:\cr 
        \code{alpha} If provided vertical lines are added, representing value 
        at risk and expected shortfall on level(s) \code{alpha} or \cr
        \code{nbins} number of supporting points, default value = 100
  }
}
\keyword{methods}