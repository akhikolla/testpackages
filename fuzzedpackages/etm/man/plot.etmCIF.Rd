\name{plot.etmCIF}
\alias{plot.etmCIF}
\title{
  Plot cumulative incidence functions
}
\description{
  Plot function for \code{etmCIF} objects. The function plots cumulative
  incidence curves, possibly with pointwise confidence intervals.
}
\usage{
\S3method{plot}{etmCIF}(x, which.cif, xlim, ylim,
    ylab = "Cumulative Incidence", xlab = "Time", col = 1, lty, lwd = 1,
    ci.type = c("none", "bars", "pointwise"), ci.fun = "cloglog",
    ci.col = col, ci.lty = 3, legend = TRUE, legend.pos, curvlab,
    legend.bty = "n", pos.ci = 27, ci.lwd = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{etmCIF} object}
  \item{which.cif}{A numeric vector indicating which CIFs should be
    plotted. When missing, only the CIF of interest is plotted
    (determined through the \code{failcode} argument in \code{\link{etmCIF}}.)}
  \item{xlim}{x-axis limits for the plot. By default, \code{c(0, max(time))}}
  \item{ylim}{y-axis limits. Default is \code{c(0, 1)}}
  \item{ylab}{Label for y-axis. Default is \code{"Cumulative Incidence"}}
  \item{xlab}{Label for x-axis. Default is "Time"}
  \item{col}{Vector describing colours used for the CIF curves. Default
    is black}
  \item{lty}{Vector of line type}
  \item{lwd}{Thickness of the lines}
  \item{ci.type}{One of \code{c("none", "bars",
      "pointwise")}. \code{none} plots no confidence interval,
    \code{bars} plots the confidence intervals in the form of a segment
    for one time point, and \code{pointwise} draws pointwise confidence
    intervals for the whole follow-up period.}
  \item{ci.fun}{Transformation used for the confidence intervals. Default is
    "clolog", and is a better choice for cumulative incidences. Other
    choices are "log" and "log-log"}
  \item{ci.col}{Colour for the pointwise confidence interval
    curves. Default is same as the CIF curves}
  \item{ci.lty}{Line type for the confidence intervals. Default is 3}
  \item{legend}{Logical. Whether to draw a legend. Default is \code{TRUE}}
  \item{legend.pos}{A vector giving the legend's position. See
    \code{\link{legend}} for further details}
  \item{curvlab}{A character or expression vector to appear in the
    legend. Default is CIF + event label}
  \item{legend.bty}{Box type for the legend. Default is none ("n")}
  \item{pos.ci}{If \code{ci.type = "bars"}, vector of integers indicating at which
    time point to put the confidence interval bars. Default is 27}
  \item{ci.lwd}{Thickness of the confidence interval segment (for
    \code{ci.type = "bars"})}
  \item{\dots}{Further graphical arguments}
}
\details{
  The function relies on \code{plot.etm} and \code{lines.etm} with more
  or less the same options. Exception is the drawing of the confidence
  intervals, for which several displays are possible.
}
\value{
  No value returned
}
\author{
  Arthur Allignol \email{arthur.allignol@gmail.com}
}
\seealso{
  \code{\link{etmCIF}}, \code{\link{plot.etm}}, \code{\link{lines.etm}}
}
\examples{
data(abortion)

cif.ab <- etmCIF(survival::Surv(entry, exit, cause != 0) ~ group, abortion,
                 etype = cause, failcode = 3)

cif.ab

plot(cif.ab, ci.type = "bars", pos.ci = 24,
     col = c(1, 2), lty = 1, curvlab = c("Control", "Exposed"))

plot(cif.ab, which = c(1, 2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
\keyword{survival}
