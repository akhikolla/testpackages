\name{aceobjects}
\alias{print.ace}
\alias{summary.ace}
\alias{print.pace}
\alias{summary.pace}

\title{
Print and Summary of \code{ace} and \code{pace} Objects
}
\description{
These functions print and summarise objects of class \code{"ace"} or \code{"pace"}.
}
\usage{
\method{print}{ace}(x, ...)
\method{summary}{ace}(object, ...)

\method{print}{pace}(x, ...)
\method{summary}{pace}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{"ace"} or \code{"pace"}.
}
  \item{object}{
An object of class \code{"ace"} or \code{"pace"}.
}
  \item{\dots}{
Arguments to be passed to and from other methods.
}
}

\value{
If the object is a result of a direct call to \code{\link{aceglm}}, \code{\link{acenlm}}, \code{\link{paceglm}}, or \code{\link{pacenlm}}, then the argument \code{criterion} will be printed, otherwise the statement \code{User-defined utility} will be printed.

Also printed are the number of repetitions (\code{"pace"} objects only), runs, factors, Phase I and II iterations of the ACE algorithm and the computational time required.

For more details on the ACE algorithm, see Overstall & Woods (2017).
}

\references{
Overstall, A.M. & Woods, D.C. (2017). Bayesian design of experiments using
approximate coordinate exchange. \emph{Technometrics}, \bold{59}, 458-470.
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}, David C. Woods, Maria Adamou & Damianos Michaelides
}

\note{For examples see \code{\link{ace}}, \code{\link{aceglm}}, and \code{\link{acenlm}}.}

\seealso{
\code{\link{ace}}, \code{\link{aceglm}}, \code{\link{acenlm}}, \code{\link{pace}}, \code{\link{paceglm}}, \code{\link{pacenlm}}.
}
