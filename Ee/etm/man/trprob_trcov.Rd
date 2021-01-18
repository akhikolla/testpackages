\name{trprob.etm}
\Rdversion{1.1}
\alias{trprob.etm}
\alias{trprob}
\alias{trcov}
\alias{trcov.etm}
\title{
  Function to extract transition probabilities and (co)variance
}
\description{
  The \code{trprob} method is used to extract transition probabilities,
  while \code{trcov} is used to obtain the (co)variance.
}
\usage{
\S3method{trprob}{etm}(x, tr.choice, timepoints, ...)
\S3method{trcov}{etm}(x, tr.choice, timepoints, ...)
}
\arguments{
  \item{x}{An object of class \code{etm}.}
  \item{tr.choice}{A character vector of the form "from to" describing
    for which transition one wishes to obtain the transition probabilities
    or covariance estimates. For \code{trprob}, \code{tr.choice} must be
    of length 1, while it can be of length 2 for \code{trcov}.}
  \item{timepoints}{Time points at which one want the estimates. When
    missing, estimates are obtained for all event times.}
  \item{\dots}{Further arguments.}
}
\value{
  A vector containing the transition probabilities or covariance
  estimates either at the time specified in \code{timepoints} or at all
  transition times.
}
\author{
  Arthur Allignol, \email{arthur.allignol@gmail.com}
}
\seealso{
  \code{\link{etm}}
}
\examples{
data(sir.cont)

# Modification for patients entering and leaving a state
# at the same date
# Change on ventilation status is considered
# to happen before end of hospital stay
sir.cont <- sir.cont[order(sir.cont$id, sir.cont$time), ]
for (i in 2:nrow(sir.cont)) {
  if (sir.cont$id[i]==sir.cont$id[i-1]) {
    if (sir.cont$time[i]==sir.cont$time[i-1]) {
      sir.cont$time[i-1] <- sir.cont$time[i-1] - 0.5
    }
  }
}

### Computation of the transition probabilities
# Possible transitions.
tra <- matrix(ncol=3,nrow=3,FALSE)
tra[1, 2:3] <- TRUE
tra[2, c(1, 3)] <- TRUE

# etm
fit.etm <- etm(sir.cont, c("0", "1", "2"), tra, "cens", 0)

## extract P_01(0, t) and variance
p01 <- trprob(fit.etm, "0 1")
var.p01 <- trcov(fit.etm, "0 1")

## covariance between P_00 and P_01
cov.00.01 <- trcov(fit.etm, c("0 0", "0 1"))

## P_01 at some time points
trprob(fit.etm, "0 1", c(0, 15, 50, 100))
}
\keyword{methods}