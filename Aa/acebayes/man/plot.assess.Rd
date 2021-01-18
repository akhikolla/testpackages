\name{plot.assess}
\alias{plot.assess}
\title{
Plot \code{assess} Objects
}
\description{
This function plots objects of class \code{"assess"}.
}
\usage{

\method{plot}{assess}(x, ...)

}
\arguments{
  \item{x}{
An object of class \code{"assess"}.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}
\value{
In the case of when \code{d1} was generated from a call to \code{(p)ace} with argument \code{deterministic = FALSE} or from a call to \code{(p)aceglm} or \code{(p)acenlm} with argument \code{method} being \code{"MC"}, then boxplots of the \code{n.assess} evaluations of the approximate expected utility for each of the designs \code{d1} and \code{d2} will be produced. Otherwise, a plot is not meaningful and a warning will be produced.
}
\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods, Maria Adamou & Damianos Michaelides
}

\seealso{
\code{\link{assess}}
}