\name{PhasePlot}
\alias{PhasePlot}
\title{Phase plot of BCPA output}
\usage{
PhasePlot(windowsweep, type = c("smooth", "flat")[1], clusterwidth = 1, ...,
  legend.where = "bottomright")
}
\arguments{
  \item{windowsweep}{\code{windowsweep} object, i.e. the
  output of the \code{\link{WindowSweep}} function.}

  \item{type}{whether to plot smooth or flat bcpa output}

  \item{clusterwidth}{for flat BCPA, this is the temporal
  range within which change points are considered to be
  within the same cluster.}

  \item{...}{additional arguments passed to the \code{plot}
  base function.}

  \item{legend.where}{where to place the tau legend (see
  \code{\link{legend}}).}
}
\description{
Behavioral phase plot of BCPA output.  Mean and standard
deviation are on the x and y axis.  Size and color of
points represent the time scale of autocorrelation
}
\examples{
if(!exists("Simp.ws"))
{
 data(Simp)
 Simp.ws <- WindowSweep(GetVT(Simp), "V*cos(Theta)", windowsize = 50, windowstep = 1, progress=TRUE)
}

PhasePlot(Simp.ws)
}
\author{
Eliezer Gurarie
}
\seealso{
\code{\link{WindowSweep}},
\code{\link{PartitionParameters}}
}

