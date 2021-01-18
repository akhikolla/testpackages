\name{plot.ace}
\alias{plot.ace}
\alias{plot.pace}
\title{
Plot \code{ace} Objects
}
\description{
This function plots objects of class \code{"ace"} or \code{"pace"} .
}
\usage{

\method{plot}{ace}(x, ...)

\method{plot}{pace}(x, ...)
}
\arguments{
  \item{x}{
An object of class \code{"ace"} or \code{"pace"}.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}
\value{
A trace plot of the current evaluations of the approximate expected utility function. Separate lines are produced for the traces from Phases I and II of the ACE algorithm. 

For objects of class \code{"pace"}, the evaluations of the approximate expected utility function are from the repetition which resulted in the terminal design (see \code{\link{pace}}, \code{\link{paceglm}}, and \code{\link{pacenlm}} for more details).

These trace plots can be used to assess convergence. See Overstall & Woods (2017) for more details.
}
\references{
Overstall, A.M. & Woods, D.C. (2017). Bayesian design of experiments using
approximate coordinate exchange. \emph{Technometrics}, \bold{59}, 458-470.
}
\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods, Maria Adamou & Damianos Michaelides
}
\note{For an example see \code{\link{ace}}.}

\seealso{
\code{\link{ace}}, \code{\link{pace}}
}