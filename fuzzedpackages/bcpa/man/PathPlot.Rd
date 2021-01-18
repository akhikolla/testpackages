\name{PathPlot}
\alias{PathPlot}
\title{Path plot of BCPA output}
\usage{
PathPlot(Data, windowsweep, type = c("smooth", "flat")[1], clusterwidth = 1,
  plotlegend = FALSE, tauwhere = "topleft", n.legend = 5,
  ncol.legend = 1, bty.legend = "n", ...)
}
\arguments{
  \item{Data}{the track data to be plotted - most
  typically, output of the \code{\link{GetVT}} function.}

  \item{windowsweep}{\code{windowsweep} object, i.e. the
  output of the \code{\link{WindowSweep}} function.}

  \item{type}{whether to plot smooth or flat bcpa output}

  \item{clusterwidth}{for flat BCPA, this is the temporal
  range within which change points are considered to be
  within the same cluster.}

  \item{plotlegend}{whether to plot a legend.}

  \item{tauwhere}{where to place the legend for the
  time-scale / auto-correlation.  Can be one of "nowhere",
  "top", "bottom", "left", "right", "topleft", "topright",
  "bottomright", "bottomleft".}

  \item{n.legend}{number of labels in legend.}

  \item{ncol.legend}{number of columns in the legend.}

  \item{bty.legend}{whether to draw a box around the
  legend.}

  \item{...}{additional arguments to pass to the
  \code{plot} base function.}
}
\description{
Plots the animal's trajectory, with segments color-coded
according to the time scale / auto-correlation of the BCPA
output, and width of segments proportional to the estimated
mean of the BCPA.
}
\examples{
if(!exists("Simp.ws"))
{
 data(Simp)
 Simp.ws <- WindowSweep(GetVT(Simp), "V*cos(Theta)", windowsize = 50, windowstep = 1, progress=TRUE)
}

PathPlot(Simp, Simp.ws, plotlegend=TRUE, n.legend=3)
PathPlot(Simp, Simp.ws, type="flat", clusterwidth=3, plotlegend=TRUE)
}
\author{
Eliezer Gurarie
}

