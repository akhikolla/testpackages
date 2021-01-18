\name{assessobjects}
\alias{print.assess}
\alias{summary.assess}

\title{
Print and Summary of \code{assess} Objects
}
\description{
These functions print and summarise objects of class \code{"assess"}.
}
\usage{
\method{print}{assess}(x, ...)
\method{summary}{assess}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{"assess"}.
}
  \item{object}{
An object of class \code{"assess"}.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}

\value{
These functions both provide a print out with the following information.

In the case of when \code{d1} was generated from a call to \code{(p)ace} with argument \code{deterministic = FALSE} or from a call to \code{(p)aceglm} or \code{(p)acenlm} with argument \code{method} being \code{"MC"}, then the mean and standard deviation of the \code{n.assess} evaluations of the approximate expected utility for each of the designs \code{d1} and \code{d2} will be printed.

Otherwise, one evaluation of the deterministic approximation to the expected utility will be printed for each of the designs from \code{d1} and \code{d2}. In the case when \code{d1} was generated as a call to \code{(p)aceglm} or \code{(p)acenlm} with argument \code{criterion} being \code{"A"}, \code{"D"} or \code{"E"}, the relative D-, E-, or A-efficiency of the two designs will be also be printed.
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods, Maria Adamou & Damianos Michaelides
}

\seealso{
\code{\link{assess}}
}
