\name{calc1DCIs}
\alias{calc1DCIs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute 1D confidence intervals
}
\description{
This computes 1D confidence intervals from an inferred likelihood surface by profile likelihood ratio methods
}
\usage{
calc1DCIs(oneDimCIvars, FONKgNames, fittedNames, CIlevel = blackbox.getOption("CIlevel"),
          nextBounds = blackbox.getOption("nextBounds"),
          NextBoundsLevel = blackbox.getOption("NextBoundsLevel"),
          boundsOutfile = "", dataString = "", cleanResu = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oneDimCIvars}{
  The names of parameters for which confidence intervals are computed
}
  \item{FONKgNames}{
  The names of \dQuote{Fitted Or Not} parameters (see Details in \code{\link{blackbox.options}} for this concept)
}
  \item{fittedNames}{
  The names of fitted parameters (see Details in \code{\link{blackbox.options}} for this concept)
}
  \item{CIlevel}{
  Level (1-coverage) of the confidence intervals. Default is 0.05.
}
  \item{nextBounds}{
    For development purposes, not documented
  }
  \item{NextBoundsLevel}{
    For development purposes, not documented
}
  \item{boundsOutfile}{
    For development purposes, not documented
}
  \item{dataString}{
  A prefix string in some outputs.
 }
  \item{cleanResu}{
  A connection, or a character string naming a file for some nicely formated output. If \code{""} (the default), print to the standard output connection.
}
}
\value{
Returns invisibly a list of profile points that met the CI level for each parameter.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Inference}
