\name{lines.etm}
\alias{lines.etm}

\title{
  Lines method for 'etm' objects
}
\description{
  Lines method for \code{etm} objects
}

\usage{
\S3method{lines}{etm}(x, tr.choice, col = 1, lty,
      conf.int = FALSE, level = 0.95, ci.fun = "linear",
      ci.col = col, ci.lty = 3, ...)
}

\arguments{
  \item{x}{An object of class \code{etm}.}
  \item{tr.choice}{character vector of the form \code{c("from to","from to")}
    specifying which transitions should be plotted. By default, all the
    direct transition probabilities are plotted}
  \item{col}{Vector of colours. Default is black.}
  \item{lty}{Vector of line type. Default is 1:number of transitions}
  \item{conf.int}{Logical specifying whether to plot confidence
  intervals. Default is FALSE.}
  \item{level}{Level of the confidence interval. Default is 0.95.}
  \item{ci.fun}{Transformation applied to the confidence intervals. It
    could be different for all transition probabilities, though if
    \code{length(ci.fun) != number of transitions}, only \code{ci.fun[1]}
    will be used. Possible choices are "linear", "log", "log-log" and
  "cloglog". Default is "linear".}
  \item{ci.col}{Colours of the confidence intervals. Default value is
    the same as \code{col}.}
  \item{ci.lty}{Line types for the confidence intervals. Default is 3.}
  \item{\dots}{Further arguments for \code{lines}.}
}

\value{
  No value returned.
}

\author{
  Arthur Allignol, \email{arthur.allignol@gmail.com}
}

\seealso{
  \code{\link{etm}}, \code{\link{plot.etm}}, \code{\link{xyplot.etm}}
}

\keyword{hplot}
\keyword{survival}
