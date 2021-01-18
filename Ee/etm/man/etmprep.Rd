\name{etmprep}
\Rdversion{1.1}
\alias{etmprep}
\title{
  Data transformation function for using etm
}
\description{
  The function transforms a data set in the wide format (i.e., one raw
  per subject) into the long format (i.e., one raw per transition, and
  possibly several raws per subjects) in a suitable way for using the
  \code{etm} function
}
\usage{
etmprep(time, status, data, tra, state.names, cens.name = NULL,
start = NULL, id = NULL, keep)
}
\arguments{
  \item{time}{A character vector giving the name of the columns
    containing the transition times or last follow-up times. The
    length of \code{time} have to be equal to the number of states, some
    elements may be NA. See Details.}
  \item{status}{A character vector giving the name of the columns
    indicating whether a state has been visited (0 if not, 1
    otherwise).}
  \item{data}{A data frame in which to look for the columns specified in
    \code{time} and \code{status}.}
  \item{tra}{A quadratic matrix of logical values describing the
    possible transitions within the multistate model. The \eqn{(i,
      j)}{(i, j)}th element of \code{tra} is TRUE if a transition from
    state \eqn{i}{i} to state \eqn{j}{j} is possible, FALSE
    otherwise. The diagonal must be set to FALSE.}
  \item{state.names}{A vector of characters giving the states names. If
    missing, state names are set to be 0:(number of states).}
  \item{cens.name}{A character string specifying how censored
    observations will be indicated in the new data set. Default is NULL,
    i.e., no censored observation.}
  \item{start}{A list containing two elements, \code{state} and
    \code{time}, giving the starting states and times for all
    individuals. Default is NULL, in which case all individuals are
    considered to start in the initial state at time 0.}
  \item{id}{A character string specifying in which column of \code{data}
    the user ids are. Default is NULL, and the ids will be \code{1:n}.}
  \item{keep}{A character vector indicating the column names of the
    covariate one might want to keep in the new data.frame.}
}
\details{
  This function only works for irreversible acyclic Markov processes. 
  Therefore, the multistate model will have initial
  states, into which no transition are possible. For these, NAs are
  allowed in \code{time} and \code{status}.
}
\value{
  The function returns a data.frame suitable for using the \code{etm}
  function. The data frame contains the following components:
  \item{id}{Individual id number}
  \item{entry}{Entry time into a state}
  \item{exit}{Exit time from a state}
  \item{from}{State from which a transition occurs}
  \item{to}{State into which a transition occurs}
  \item{\dots}{Further columns specified in \code{keep}}
}
\author{
  Arthur Allignol, \email{arthur.allignol@gmail.com}
}
\seealso{
  \code{\link{etm}}
}
\examples{
### creation of fake data in the wild format, following an illness-death model
## transition times
tdisease <- c(3, 4, 3, 6, 8, 9)
tdeath <- c(6, 9, 8, 6, 8, 9)

## transition status
stat.disease <- c(1, 1, 1, 0, 0, 0)
stat.death <- c(1, 1, 1, 1, 1, 0)

## a covariate that we want to keep in the new data
cova <- rbinom(6, 1, 0.5)

dat <- data.frame(tdisease, tdeath,
                  stat.disease, stat.death,
                  cova, stringsAsFactors = TRUE)

## Possible transitions
tra <- matrix(FALSE, 3, 3)
tra[1, 2:3] <- TRUE
tra[2, 3] <- TRUE

## data preparation
newdat <- etmprep(c(NA, "tdisease", "tdeath"),
                  c(NA, "stat.disease", "stat.death"),
                  data = dat, tra = tra, cens.name = "cens")
}
\keyword{datagen}
\keyword{manip}
