\name{plot.bcpa}
\alias{plot.bcpa}
\title{Plotting method for BCPA output}
\usage{
\method{plot}{bcpa}(x, type = c("smooth", "flat")[1], threshold = 3,
  clusterwidth = 1, col.cp = rgb(0.5, 0, 0.5, 0.5), pt.cex = 0.5,
  legend = TRUE, rho.where = "topleft", mu.where = "nowhere",
  col.sd = "red", col.mean = "black", ...)
}
\arguments{
  \item{x}{a \code{windowsweep} object, i.e. the output of
  the \code{\link{WindowSweep}} function.}

  \item{type}{whether to plot smooth or flat bcpa output}

  \item{threshold}{for smooth BCPA, this is the minimum
  number of windows that must have identified a given
  changepoint to be illustrated.}

  \item{clusterwidth}{for flat BCPA, this is the temporal
  range within which change points are considered to be
  within the same cluster.}

  \item{col.cp,col.mean,col.sd}{color of the vertical
  change points, mean estimate, and prediction interval (mu
  +- sigma), respectively.}

  \item{pt.cex}{expansion coefficient for point sizes.}

  \item{legend}{logical - whether to draw a legend or not.}

  \item{rho.where}{where to place the legend for the
  time-scale / auto-correlation.  Can be one of "nowhere",
  "top", "bottom", "left", "right", "topleft", "topright",
  "bottomright", "bottomleft"}

  \item{mu.where}{where (and whether) to place the legend
  box for the mean - same options as for \code{rho.where}}

  \item{...}{other arguments to pass to the \code{plot}
  base function.}
}
\description{
Plotting method for the output of a BCPA analysis with
vertical break points, superimposed estimates of the
partitioned mean and variance estimates and color-coded
autocorrelation estimates.
}
\examples{
if(!exists("Simp.ws"))
{
 data(Simp)
 Simp.ws <- WindowSweep(GetVT(Simp), "V*cos(Theta)", windowsize = 50, windowstep = 1, progress=TRUE)
}

plot(Simp.ws)
# this actually provides basically the exact original changepoints
plot(Simp.ws, threshold=7)
# here's a flat analysis
plot(Simp.ws, type="flat", clusterwidth=3, legend=FALSE)
}
\author{
Eliezer Gurarie
}
\seealso{
Plots output of the \code{\link{WindowSweep}} function.
}

