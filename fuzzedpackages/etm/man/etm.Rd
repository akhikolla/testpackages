\name{etm}
\alias{etm}
\alias{etm.data.frame}

\title{Computation of the empirical transition matrix}
\description{
  This function computes the empirical transition matrix, also called
  Aalen-Johansen estimator, of the transition probability matrix of any
  multistate model. The covariance matrix is also computed.
}
\usage{
\S3method{etm}{data.frame}(data, state.names, tra, cens.name, s, t = "last",
    covariance = TRUE, delta.na = TRUE, modif = FALSE,
    c = 1, alpha = NULL, strata, ...)
}
\arguments{
  \item{data}{data.frame of the form data.frame(id,from,to,time)
    or (id,from,to,entry,exit)
    \describe{
    \item{id:}{patient id}
    \item{from:}{the state from where the transition occurs}
    \item{to:}{the state to which a transition occurs}
    \item{time:}{time when a transition occurs}
    \item{entry:}{entry time in a state}
    \item{exit:}{exit time from a state}
    }
  This data.frame is transition-oriented, \emph{i.e.} it contains one
  row per transition, and possibly several rows per patient. Specifying
  an entry and exit time permits to take into account left-truncation. }
  \item{state.names}{A vector of characters giving the states names.}
  \item{tra}{A quadratic matrix of logical values describing the possible
    transitions within the multistate model. }
  \item{cens.name}{ A character giving the code for censored
    observations in the column 'to' of \code{data}. If there is no
    censored observations in your data, put 'NULL'.}
  \item{s}{Starting value for computing the transition probabilities.}
  \item{t}{Ending value. Default is "last", meaning that the transition
    probabilities are computed over \eqn{(s, t]}{(s, t]}, \eqn{t}{t}
    being the last time in the data set.}
  \item{covariance}{Logical. Decide whether or not computing the
    covariance matrix. May be useful for, say, simulations, as the variance
    computation is a bit long. Default is TRUE.}
  \item{delta.na}{Logical. Whether to export the array containing the
    increments of the Nelson-Aalen estimator. Default is \code{TRUE}.}
  \item{modif}{Logical. Whether to apply the modification of Lai and
    Ying for small risk sets}
  \item{c}{Constant for the Lai and Ying modification. Either \code{c}
    contains only one value that will be used for all the states,
    otherwise \code{c} should be the same length as
    \code{state.names}.}
  \item{alpha}{Constant for the Lai and Ying modification. If NULL (the
    default) then only \code{c} is used and the Lai and Ying
    modification discards the event times for which \eqn{Y(t) \geq
      t}{Y(t) >= t}. Otherwise \eqn{cn^\alpha}{cn^alpha} is used. It is
    recommanded to let \code{alpha} equal NULL for multistate models.}
  \item{strata}{Character vector giving variables on which to stratify
    the analysis.}
  \item{...}{Not used}
}
\details{
  Data are considered to arise from a time-inhomogeneous Markovian
  multistate model with finite state space, and possibly subject to
  independent right-censoring and left-truncation.

  The matrix of the transition probabilities is estimated by the
  Aalen-Johansen estimator / empirical transition matrix (Andersen et
  al., 1993), which is the product integral over the time period
  \eqn{(s, t]}{(s, t]} of I + the matrix of the increments of the
  Nelson-Aalen estimates of the cumulative transition hazards. The
  \eqn{(i, j)-th}{(i, j)-th} entry of the empirical transition matrix
  estimates the transition probability of being in state \eqn{j}{j} at
  time \eqn{t}{t} given that one has been in state j at time \eqn{s}{s}.

  The covariance matrix is computed using the recursion formula (4.4.19)
  in Anderson et al. (1993, p. 295). This estimator of the covariance
  matrix is an estimator of the Greenwood type.

  If the multistate model is not Markov, but censorship is entirely
  random, the Aalen-Johansen estimator still consistently estimates the
  state occupation probabilities of being in state \eqn{i}{i} at time
  \eqn{t}{t} (Datta & Satten, 2001; Glidden, 2002)

  Recent versions of R have changed the \code{data.frame} function,
  where the default for the \code{stringsAsFactors}
  argument from \code{TRUE} to \code{FALSE}. \code{etm} currently
  depends on the states being factors, so that the user should use
  \code{data.frame(..., stringsAsFactors=TRUE)}.
  
  }
\value{
  \item{est}{Transition probability estimates. This is a 3 dimension
    array with the first dimension being the state from where transitions
    occur, the second the state to which transitions occur, and the
    last one being the event times.}
  \item{cov}{Estimated covariance matrix. Each cell of the matrix gives
    the covariance between the transition probabilities given by the
    rownames and the colnames, respectively.}
  \item{time}{Event times at which the transition probabilities are
  computed. That is all the observed times between \eqn{(s, t]}{(s, t]}.}
  \item{s}{Start of the time interval.}
  \item{t}{End of the time interval.}
  \item{trans}{A \code{data.frame} giving the possible transitions.}
  \item{state.names}{A vector of character giving the state names.}
  \item{cens.name}{How the censored observation are coded in the data
    set.}
  \item{n.risk}{Matrix indicating the number of individuals at risk just
    before an event}
  \item{n.event}{Array containing the number of transitions at each
    times}
  \item{delta.na}{A 3d array containing the increments of the
    Nelson-Aalen estimator.}
  \item{ind.n.risk}{When \code{modif} is true, risk set size for which
    the indicator function is 1}

  If the analysis is stratified, a list of \code{etm} objects is
  returned.
}
\references{
  Beyersmann J, Allignol A, Schumacher M: Competing Risks and Multistate
  Models with R (Use R!), Springer Verlag, 2012 (Use R!)
  
  Allignol, A., Schumacher, M. and Beyersmann, J. (2011).
  Empirical Transition Matrix of Multi-State Models: The etm Package.
  \emph{Journal of Statistical Software}, 38.

  Andersen, P.K., Borgan, O., Gill, R.D. and Keiding,
  N. (1993). \emph{Statistical models based on counting
    processes}. Springer Series in Statistics. New York, NY: Springer.
  
  Aalen, O. and Johansen, S. (1978). An empirical transition matrix for
  non-homogeneous Markov chains based on censored
  observations. \emph{Scandinavian Journal of Statistics}, 5: 141-150.

  Gill, R.D. and Johansen, S. (1990). A survey of product-integration
  with a view towards application in survival analysis. \emph{Annals of
    statistics}, 18(4): 1501-1555.

  Datta, S. and Satten G.A. (2001). Validity of the Aalen-Johansen
  estimators of stage occupation probabilities and Nelson-Aalen
  estimators of integrated transition hazards for non-Markov
  models. \emph{Statistics and Probability Letters}, 55(4): 403-411.

  Glidden, D. (2002). Robust inference for event probabilities with
  non-Markov data. \emph{Biometrics}, 58: 361-368.
}
\author{Arthur Allignol, \email{arthur.allignol@gmail.com}}

\note{Transitions into a same state, mathematically superfluous, are not
  allowed. If transitions into the same state are detected in the data,
  the function will stop. Equally, \code{diag(tra)} must be set to
  FALSE, see the example below.}

\seealso{\code{\link{print.etm}}, \code{\link{summary.etm}}, \code{\link{sir.cont}},
  \code{\link{xyplot.etm}}}

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
tr.prob <- etm(sir.cont, c("0", "1", "2"), tra, "cens", 1)

tr.prob
summary(tr.prob)

# plotting
if (require("lattice")) {
xyplot(tr.prob, tr.choice=c("0 0", "1 1", "0 1", "0 2", "1 0", "1 2"),
       layout=c(2, 3), strip=strip.custom(bg="white",
         factor.levels=
     c("0 to 0", "1 to 1", "0 to 1", "0 to 2", "1 to 0", "1 to 2")))
}

### example with left-truncation

data(abortion)

# Data set modification in order to be used by etm
names(abortion) <- c("id", "entry", "exit", "from", "to")
abortion$to <- abortion$to + 1

## computation of the matrix giving the possible transitions
tra <- matrix(FALSE, nrow = 5, ncol = 5)
tra[1:2, 3:5] <- TRUE

## etm
fit <- etm(abortion, as.character(0:4), tra, NULL, s = 0)

## plot
xyplot(fit, tr.choice = c("0 0", "1 1", "0 4", "1 4"),
       ci.fun = c("log-log", "log-log", "cloglog", "cloglog"),
       strip = strip.custom(factor.levels = c("P(T > t) -- control",
                                              "P(T > t) -- exposed",
                                 "CIF spontaneous abortion -- control",
                                 "CIF spontaneous abortion --
exposed")))
}

\keyword{survival}
