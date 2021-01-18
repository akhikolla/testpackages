\name{tra}
\alias{tra}
\alias{tra_ill}
\alias{tra_ill_comp}
\alias{tra_comp}
\alias{tra_surv}

\title{
  Matrix of possible transitions
}
\description{
  Miscellaneous functions that compute the matrix of possible
  transitions used as argument in the \code{etm} function.
}
\usage{
tra_ill(state.names = c("0", "1", "2"))
tra_ill_comp(nComp = 2,
     state.names = as.character(seq(0, nComp + 1, 1)))
tra_comp(nComp = 2,
         state.names = as.character(seq(0, nComp)))
tra_surv(state.names = c("0", "1"))
}

\arguments{
  \item{state.names}{A vector of characters giving the states names}
  \item{nComp}{For the competing risks models, the number of competing
  events}
}
\details{
  These functions compute the matrix of possible transitions that is
  used as argument in, e.g., the \code{etm} function. \code{tra_surv} is
  for the usual survival model, \code{tra_comp} for the competing risks
  model, \code{tra_ill} for the illness-death model and
  \code{tra_ill_comp} for the illness-death model with competing
  terminal events. By default, state names are from 0 to \dots
}
\value{
  A quadratic matrix with \code{TRUE} if a transition is possible,
  \code{FALSE} otherwise.
}

\author{
  Arthur Allignol \email{arthur.allignol@gmail.com}
}

\seealso{
  \code{\link{etm}}
}
\examples{
tra_ill()

## competing risks model with 4 competing events non-default state names
tra_comp(4, state.names = c("healthy", "Cardiac problems", "Cancer",
                           "Rhenal failure", "Other"))
}

\keyword{survival}
\keyword{miscellaneous}
