\name{GetModels}
\alias{GetModels}
\title{Model selection at a known breakpoint}
\usage{
GetModels(x, t, breakpoint, K = 2, tau = TRUE)
}
\arguments{
  \item{x}{vector of time series values.}

  \item{t}{vector of times of measurements associated with
  x.}

  \item{breakpoint}{breakpoint to test (in terms of the
  INDEX within "t" and "x", not actual time value).}

  \item{K}{sensitivity parameter.  Standard definition of
  BIC, K = 2.  Smaller values of K mean less sensitive
  selection, i.e. higher likelihood of selecting null (or
  simpler) models.}

  \item{tau}{whether or not to estimate time scale
  \eqn{\tau} (preferred) or autocorrelation \eqn{\rho}.}
}
\value{
Returns a names matrix with 8 rows (one for each model) and
columns: \code{Model}, \code{LL, bic, mu1, s1, rho1, mu2,
s2, rho2}.  Fairly self-explanatory.  Note that the
\code{rho} columns include the \code{tau} values, if
\code{tau} is TRUE (as it usually should be).
}
\description{
Returns all parameter estimates and log-likelihoods for all
possible models at a selected breakpoint.  These are:
\itemize{ \item{M0} - all parameters equal \item{M1} -
\eqn{\mu_1 != \mu_2} \item{M2} - \eqn{\sigma_1 != \sigma_2}
\item{M3} - \eqn{\tau_1 != \tau_2} \item{M4} - \eqn{\mu_1
!= \mu_2} AND \eqn{\sigma_1 != \sigma_2} \item{M5} -
\eqn{\mu_1 != \mu_2} AND \eqn{\tau_1 != \tau_2} \item{M6} -
\eqn{\sigma_1 != \sigma_2} AND \eqn{\tau_1 != \tau_2}
\item{M7} - all parameters unequal}
}
\author{
Eliezer Gurarie
}
\seealso{
Used directly within \code{\link{WindowSweep}}.  Relies
heavily on \code{\link{GetRho}}.
}

