\name{WindowSweep}
\alias{WindowSweep}
\title{Perform window sweep for BCPA}
\usage{
WindowSweep(data, variable, windowsize = 50, windowstep = 1, K = 2,
  tau = TRUE, range = 0.6, progress = TRUE, plotme = FALSE, ...)
}
\arguments{
  \item{data}{the data to be analyzed.  Most typically,
  output of the \code{\link{GetVT}} function containing
  step lengths, absolute and turning angles, etc.}

  \item{variable}{a character string representing the
  response to apply the BCPA to.  For example
  \code{"V*cos(theta)"} for persistence velocity, or
  \code{"log(V)"} for log of velocity.}

  \item{windowsize}{integer size of the analysis window as
  a number of data points (not time units).  Should
  probably be no smaller than 20.}

  \item{windowstep}{integer step size of analysis. Values
  greater than 1 speed the analysis up.}

  \item{K}{senstivity parameter for the adjusted BIC.
  Smaller values make for a less sensitive model selection,
  i.e. more likely that the null model of no significant
  changes will be selected.}

  \item{tau}{a logical indicating whether the
  autocorrelation "rho" or the characteristic time "tau"
  should be estimated.}

  \item{range}{a number between 0 and 1 that determines the
  extent of each window that is scanned for changepoints.
  I.e., if the window is 100 datapoints long, at the
  default \code{range=0.6}, changepoints will be scanned
  between 20 and 80.}

  \item{progress}{logical - whether or not to output a
  progress bar as the analysis is being performed.}

  \item{plotme}{logical - whether or not to plot the
  analysis as it is happening.  This slows the analysis
  considerably, especially in non-dynamic graphic
  environments like RStudio.}

  \item{...}{additional parameters to be passed to the
  \code{\link{PartitionParameters}} function.}
}
\value{
an object of class \code{windowsweep}, which is a list
containing: \item{ws}{a data frame containing the change
point, selected model, and parameter estimates}
\item{x}{the response variable} \item{t}{the time variable
- in the units specified in the data object}
\item{t.POSIX}{the time variable as a POSIX objects
(contianing Y-M-D H:S)} \item{windowsize}{the window size}
\item{windowstep}{the window step size}
}
\description{
This is the workhorse function of the BCPA.  It performs a
sweep of the time series, searching for most significant
change points and identifying the parsimonious model
according to an adjusted BIC.
}
\examples{
data(Simp)
plot(Simp)
Simp.VT <- GetVT(Simp)
Simp.ws <- WindowSweep(Simp.VT, "V*cos(Theta)", windowsize = 50, windowstep = 1, progress=TRUE)
plot(Simp.ws, threshold=7)
plot(Simp.ws, type="flat", clusterwidth=3)
PathPlot(Simp, Simp.ws)
PathPlot(Simp, Simp.ws, type="flat")
DiagPlot(Simp.ws)
}
\author{
Eliezer Gurarie
}
\seealso{
for internal functions: \code{\link{GetModels}},
\code{\link{GetBestBreak}}, \code{\link{GetDoubleL}}; for
summarizing output: \code{\link{ChangePointSummary}}; for
plotting output: \code{\link{plot.bcpa}}
}

