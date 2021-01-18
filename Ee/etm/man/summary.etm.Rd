\name{summary.etm}
\alias{summary.etm}
\alias{print.summary.etm}
\title{Summary methods for an 'etm' object}
\description{
  Summary method for objects of class \code{etm}
}
\usage{
\S3method{summary}{etm}(object, tr.choice,
        ci.fun = "linear", level = 0.95, times, ...)
\S3method{print}{summary.etm}(x, ...)
}
\arguments{
  \item{object}{An object of class \code{etm}.}
  \item{tr.choice}{Character vector of the form 'c("from to","from
    to")' specifying which transitions should be summarized. Default to all
    the transition probabilities}
  \item{ci.fun}{A character vector specifying the transformation to be
    applied to the pointwise confidence intervals. It could be different
    for each transition probability, though if \code{length(ci.fun) !=
      number of transitions}, only \code{ci.fun[1]} will be used. The
    function displays the transition probabilities in the following
    order: first the direct transitions in alphabetical order, e.g., 0 to
    1, 0 to 2, 1 to 2, ..., then the state occupation probabilities in
    alphabetical order, e.g., 0 to 0, 1 to 1, ... The
    possible transformations are "linear", "log", "log-log" and
    "cloglog". Default is "linear".}
  \item{level}{Level of the two-sided confidence intervals. Default is
    0.95.}
  \item{x}{A \code{summary.cpf} object}
  \item{times}{Time points for which estimates should be
    returned. Default to all transition times.}
  \item{\dots}{Further arguments}
}
\value{
  A list of data.frames giving the transition probability and stage
  occupation probability estimates. List items are named after the
  possible transition.
  \item{P}{Transition probability estimates}
  \item{var}{Variance estimates}
  \item{lower}{Lower confidence limit}
  \item{upper}{Upper confidence limit}
  \item{time}{Transition times}
  \item{n.risk}{Number of individuals at risk of experiencing a transition
    just before time \eqn{t}{t}}
  \item{n.event}{Number of events at time \eqn{t}{t}}
}
\author{Arthur Allignol \email{arthur.allignol@gmail.com}}
\seealso{\code{\link{etm}}}
\keyword{methods}
\keyword{print}
