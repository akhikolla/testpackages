\name{GetBestBreak}
\alias{GetBestBreak}
\title{Find most likely change point in irregular time series}
\usage{
GetBestBreak(x, t, range = 0.6, ...)
}
\arguments{
  \item{x}{vector of time series values.}

  \item{t}{vector of times of measurements associated with
  x.}

  \item{range}{tange of possible breaks. Default (0.6) runs
  approximately from 1/5th to 4/5ths of the total length of
  the time series.}

  \item{...}{additional parameters to pass to
  \code{\link{GetDoubleL}} function.}
}
\value{
returns a single row (vector) with elements:
\code{breaks},\code{tbreaks},\code{mu1},\code{sigma1},\code{rho1},\code{LL1},\code{mu2},\code{sigma2},\code{rho2},\code{LL2},\code{LL}.
The breakpoint is calculated for a range of possible values
of width \code{range*l} (where \code{l} is the length of
the time series). The output of this function feeds
\code{\link{WindowSweep}}.
}
\description{
Finds the single best change point according to the
likelihood function.  Used internally within
\code{\link{WindowSweep}}.
}
\examples{
# An example with a single break:
x <- c(arima.sim(list(ar = 0.9), 20) + 10, arima.sim(list(ar = 0.1), 20))
t <- 1:length(x)
plot(t,x, type="l")
(bb <- GetBestBreak(x,t, tau=FALSE))
abline(v = bb[2], col=2)
}
\author{
Eliezer Gurarie
}
\seealso{
\code{\link{WindowSweep}} which uses it, and
\code{\link{GetDoubleL}} for the likelihood estimation.
}

