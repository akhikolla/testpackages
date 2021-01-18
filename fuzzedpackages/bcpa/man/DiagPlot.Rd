\name{DiagPlot}
\alias{DiagPlot}
\title{Diagnostic plot for BCPA}
\usage{
DiagPlot(windowsweep, type = c("smooth", "flat")[1], plotme = TRUE,
  values = FALSE, ...)
}
\arguments{
  \item{windowsweep}{a \code{windowsweep} object, i.e. the
  output of the \code{\link{WindowSweep}} analysis.}

  \item{type}{whether to diagnose the model fitted for a
  smooth or flat BCPA.}

  \item{plotme}{logical - whether or not to plot the
  diagnostics}

  \item{values}{logical - whether or not to return the
  values of the standardized residuals.}

  \item{...}{additional arguments to pass to the
  \code{\link{PartitionParameters}} function.}
}
\value{
If \code{values} is TRUE, returnn the values of the
standardized residuals.
}
\description{
Draws diagnostic plots for BCPA analysis.  Specifically: a
qqplot, a histogram (with a N(0,1) density curve), and an
acf of the standardized residuals of an analysis.
}
\examples{
data(Simp)
if(!exists("Simp.VT"))
 Simp.VT <- GetVT(Simp)
if(!exists("Simp.ws"))
 Simp.ws <- WindowSweep(Simp.VT, "V*cos(Theta)", windowsize = 50, windowstep = 1, progress=TRUE)
DiagPlot(Simp.ws)
DiagPlot(Simp.ws, type="flat")
# The Simp's diagnostic plots are excellent.
}
\author{
Eliezer Gurarie
}
\seealso{
\code{\link{PartitionParameters}}
}

