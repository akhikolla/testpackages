\name{prepare.los.data}
\alias{prepare.los.data}
\title{Prepare the data for clos}
\description{Prepare data to be passed to clos() in package etm.}
\usage{prepare.los.data(x) }
\arguments{
  \item{x}{data.frame of the form data.frame(id, j.01, j.02, j.03, j.12, j.13, cens):
    \describe{
    \item{id:}{id (patient id, admision id)}      
    \item{j.01:}{observed time for jump from 0 to 1}
    \item{j.02:}{observed time for jump from 0 to 2}
    \item{j.03:}{observed time for jump from 0 to 3}
    \item{j.12:}{observed time for jump from 1 to 2}
    \item{j.13:}{observed time for jump from 1 to 3}
    \item{cens:}{censoring time (either in initial or intermediate state)}
  }
 }
}
\value{ a data.frame of the form data.frame(id, from, to, time, oid):
  \item{id:}{ id (patient id, admision id)}
  \item{from:}{ the state from where a transition occurs}
  \item{to:}{ the state to which a transition occurs}
  \item{time:}{ time of the transition}
  \item{oid:}{ the observation id}
}

\author{ Matthias Wangler}
\seealso{ \code{\link[etm]{clos}}}
\examples{
data(los.data)
my.observ <- prepare.los.data(x=los.data)
}

\keyword{datasets}
\keyword{manip}
