\name{ChangePointSummary}
\alias{ChangePointSummary}
\title{Obtain summary of BCPA analysis}
\usage{
ChangePointSummary(windowsweep, clusterwidth = 1, tau = TRUE)
}
\arguments{
  \item{windowsweep}{a \code{windowsweep} object, i.e. the
  output of the \code{\link{WindowSweep}} function.}

  \item{clusterwidth}{the temporal range within which
  change points are considered to be within the same
  cluster.  Corresponds to the bandwidth of the density of
  the break distribution.}

  \item{tau}{logical, whether to estimate the
  characteristic time tau (preferred) or not.  If FALSE,
  the autocorrelation parameter rho is calculated.}
}
\value{
a list containing two elements: \item{breaks}{a data frame
containing columns: \code{middle} - each change point,
\code{size} - the number of windows that selected the
change point, \code{modelmode} - the most frequently
selected of the seven possible models (M0 - the null model
- is excluded), and \code{middle.POSIX} - the mid-point as
a POSIX time object.} \item{phases}{a data frame containing
columns \code{mu.hat, s.hat, rho.hat} - the estimated mean,
standard deviation, and time scale (or auto-correlation)
within each phase (i.e. period between change points),
\code{t0} - the beginning of the phase, \code{t1} - the end
of the phase, and \code{interval} - the total duration of
the phase.}
}
\description{
Produces a summary of change points for a "flat" analysis,
identifying phases (periods between change points) with
estimated parameters, clustering neighboring ones according
to a kernel density of the windowsweep breaks.
}
\examples{
if(!exists("Simp.VT")){
 data(Simp)
 Simp.VT <- GetVT(Simp)}
if(!exists("Simp.ws"))
 Simp.ws <- WindowSweep(Simp.VT, "V*cos(Theta)", windowsize = 50, windowstep = 1, progress=TRUE)
# too many change points:
ChangePointSummary(Simp.ws)
# about the right number of change points:
ChangePointSummary(Simp.ws, clusterwidth=3)
}
\author{
Eliezer Gurarie
}
\seealso{
\code{\link{WindowSweep}}
}

