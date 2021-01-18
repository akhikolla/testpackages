\name{GetVT}
\alias{GetVT}
\title{Obtain VT table from Track}
\usage{
GetVT(Data, units = "hour", skiplast = TRUE)
}
\arguments{
  \item{Data}{a track to analyze.  Must contain columns: X,
  Y and Time (as a POSIX object).  The \code{track} class
  is a robust entry.}

  \item{units}{the time units for the analysis; one of
  \code{sec, min, hour, day}.}

  \item{skiplast}{filters away last step.}
}
\value{
a data frame with the following columns: \item{Z.start,
Z.end}{the start and end locations (as complex
coordinates)} \item{S}{the step length} \item{Phi,
Theta}{absolute and turning angle, respectively}
\item{T.start, T.end}{start and time of steps (numeric - in
given units)} \item{T.mid}{temporal midpoint of the step }
\item{dT}{duration of the step} \item{V}{approximate speed
(S/dT)} \item{T.POSIX}{the temporal midpoint of the step as
a POSIX objects.}
}
\description{
The VT table containes speeds, steplengths, orientations
and other summaries derived from a track. The output of
this function is (typically) meant to feed the
\code{\link{WindowSweep}} function.
}
\examples{
data(Simp)
plot(Simp)
Simp.VT <- GetVT(Simp)
head(Simp.VT)
# Distribution of estimated speeds
hist(Simp.VT$V, col="grey", breaks=20)
# Distribution of turning angles
require(circular)
rose.diag(Simp.VT$Theta, bins=24)
}
\author{
Eliezer Gurarie
}
\seealso{
\code{\link{WindowSweep}}
}

