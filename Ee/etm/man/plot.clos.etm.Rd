\name{plot.clos.etm}
\Rdversion{1.1}
\alias{plot.clos.etm}
\title{
  Plot method for 'clos.etm' objects
}
\description{
  Plot method for objects of class \code{clos.etm}.
}
\usage{
\S3method{plot}{clos.etm}(x, xlab = "Time", ylab.e = "Expected LOS",
ylab.w = "Weights", xlim, ylim.e, ylim.w, col.e = c(1, 2), col.w = 1,
lty.e = c(1, 1), lty.w = 1, legend = TRUE, legend.pos, curvlab,
legend.bty = "n", ...)
}
\arguments{
  \item{x}{An object of class \code{clos.etm}}
  \item{xlab}{Label for the x-axis}
  \item{ylab.e}{Label for the y-axis in the plot of the expected LOS}
  \item{ylab.w}{Label for the y-axis in the plot of the weights}
  \item{xlim}{Limits of x-axis for the plots}
  \item{ylim.e}{Limits of the y-axis for the expected LOS plot}
  \item{ylim.w}{Limits of the y-axis for the weights plot}
  \item{col.e}{Vector of colours for the plot of expected LOS}
  \item{col.w}{Vector of colours for the plot of the weights}
  \item{lty.e}{Vector of line type for the plot of expected LOS}
  \item{lty.w}{Vector of line type for the plot of the weights}
  \item{legend}{Logical. Whether to draw a legend for the plot of
    expected LOS}
  \item{legend.pos}{A vector giving the legend's position. See
    \code{\link{legend}} for details}
  \item{curvlab}{Character or expression vector to appear in the
    legend. Default is \code{c("Intermediate event by time t", "No
      intermediate event by time t")}}
  \item{legend.bty}{Box type for the legend}
  \item{\dots}{Further arguments for plot}
}
\details{
  Two graphs are drawn. The lower graph displays the expected LOS for
  patients who have experienced the intermediate event and for those who
  have not. The upper graph displays the weights used to compute
  the weighted average.
}
\value{
  No value returned
}
\author{
  Arthur Allignol \email{arthur.allignol@gmail.com}, Matthias Wangler
}
\seealso{
  \code{\link{clos}}
}
\keyword{hplot}