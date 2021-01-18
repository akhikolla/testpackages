\name{bcpa-package}
\alias{bcpa-package}
\alias{bcpa}
\docType{package}
\title{Behavioral Change Point Analysis}
\description{

A collection of functions that allows one to perform the behavioral change point analysis (BCPA) as described by Gurarie et al. (2009, Ecology Letters, 12: 395-408).  The key features are estimation of discrete changes in time-series data, notable linear and turning components of gappy velocity times series extracted from movement data.  

There is a fairly detailed vignette accesible by entering \code{vignette("bcpa")}.  Alternatively, the key analysis function is \code{\link{WindowSweep}}, and reading its documentation is a good way to start using this package. This function uses a suite of functions that might also be useful for more narrow analysis, listed hierarchically (from bottom-up) below:
	\tabular{ll}{
    \code{\link{GetRho}} \tab maximizes the likelihood to estimate autocorrelation rho or characteristic time-scale tau.\cr
		\code{\link{GetDoubleL}} \tab estimates the paramters and returns the log-likelhood at either side of a given break\cr
    \code{\link{GetBestBreak}} \tab finds the single best change point according to the likelihood returned by \code{\link{GetDoubleL}} \cr    
    \code{\link{GetModels}} \tab uses a (modified) BIC model selection for all combinations from M0 (\eqn{\mu_1 = \mu_2}, \eqn{\sigma_1 = \sigma_2}, \eqn{\rho_1 = \rho_2}) to M7 (\eqn{\mu_1 \neq \mu_2}, \eqn{\sigma_1 \neq \sigma_2}, \eqn{\rho_1 \neq \rho_2}) to characterize the "Best Break"\cr
    \code{\link{WindowSweep}} \tab sweeps a longer time series with the Best Break / Model Selection analysis, identifying most likely break points and BIC selected models across the time series. \cr
      }
Summary, diagnostic, and plotting functions are:
  \tabular{ll}{
\code{\link{PartitionParameters}} \tab outputs the estimated parameters of a bcpa.\cr
\code{\link{ChangePointSummary}} \tab provides a summary table of the chage points.\cr
\code{\link{plot.bcpa}} \tab a plotting method for visualizing the time series with vertical lines as change points.\cr
\code{\link{PathPlot}} \tab a method for drawing a color-coded track of the analysis.\cr
\code{\link{DiagPlot}} \tab diagnostic plots for BCPA.\cr
	}
  
A few preprocessing functions available:
  \tabular{ll}{
\code{\link{plot.track}} \tab method for plotting a generic "track" object.\cr
\code{\link{GetVT}} \tab returns step-lenghts, absolute and turning angles from track data.\cr
  }

}
\details{
\tabular{ll}{
Package: \tab bcpa\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-10-23\cr
License: \tab Unlimited\cr
LazyLoad: \tab yes\cr
}
This is a suite of functions needed to perform a complete behavioral change point analysis.  
}
\author{Eliezer Gurarie <eliezg@uw.edu>}
\references{
	Gurarie, E., R. Andrews and K. Laidre. 2009. A novel method for identifying behavioural changes in animal movement data. Ecol. Lett. 12: 395-408.
}

\examples{
# Running through a complete analysis here:
## loading the data
  data(Simp)
## plotting the track (using the plot.track method)
  plot(Simp)
## Obtaining the movement summary table (with turning angles and step lengths)
  Simp.VT <- GetVT(Simp)
## Applying the analysis
  Simp.ws <- WindowSweep(Simp.VT, "V*cos(Theta)", windowsize = 50, windowstep = 1, progress=TRUE)
## plotting outpots
  plot(Simp.ws, threshold=7)
  plot(Simp.ws, type="flat", clusterwidth=3)
  PathPlot(Simp, Simp.ws)
  PathPlot(Simp, Simp.ws, type="flat")
## Diagnostic of assumptions
  DiagPlot(Simp.ws)
}


\keyword{bcpa}
