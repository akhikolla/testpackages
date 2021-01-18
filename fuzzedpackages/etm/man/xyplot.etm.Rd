\name{xyplot.etm}
\alias{xyplot.etm}
\title{xyplot method for object of class 'etm'}
\description{
  xyplot function for objects of class \code{etm}. Estimates of the 
  transition probabilities are plotted as a function of time for all the
  transitions specified by the user.
}
\usage{
\S3method{xyplot}{etm}(x, data = NULL, tr.choice, col = c(1, 1, 1),
       lty = c(1, 3, 3),  xlab = "Time",
       ylab = "Transition probability",
       conf.int = TRUE, ci.fun = "linear", level = 0.95, ...)
}
\arguments{
  \item{x}{An object of class \code{etm}.}
  \item{data}{\emph{Useless}.}
  \item{tr.choice}{A character vector of the form c("from to", "from
    to", ...) specifying the transition probabilities to be plotted. By
    default, all the direct transition probabilities are displayed.}
  \item{col}{Vector of colours for the curves.}
  \item{lty}{Vector of line types.}
  \item{xlab}{x-axis label. Default is "Time".}
  \item{ylab}{y-axis label. Default is "Estimated transition
    probability".}
  \item{conf.int}{Logical. Whether to draw pointwise confidence
    intervals. Default is TRUE.}
  \item{ci.fun}{A character vector specifying the transformation to be
    applied to the pointwise confidence intervals. It could be different
    for each transition probability, though if \code{length(ci.fun) !=
      length(tr.choice)}, only \code{ci.fun[1]} will be used. The
    possible transformations are "linear", "log", "log-log" and
    "cloglog". Default is "linear".}
  \item{level}{Level of the two-sided confidence intervals. Default is
    0.95.}
  \item{\dots}{Further arguments for \code{xyplot}.}
}
\value{
  An object of class \code{trellis}.
}
\author{Arthur Allignol, \email{arthur.allignol@gmail.com}}
\seealso{\code{\link{etm}}, \code{\link[lattice]{xyplot}}}
\keyword{hplot}
