% vi:et:ts=2
\name{Ropj-package}
\alias{Ropj-package}
\alias{Ropj}
\docType{package}
\title{\packageTitle{Ropj}}
\description{\packageDescription{Ropj}}
\details{
  The DESCRIPTION file: \packageDESCRIPTION{Ropj}
  \packageIndices{Ropj}

  This package exports a single function, \command{read.opj}, that tries
  to read an OPJ file and return a list of objects that it consists of.

  No write support seems to be planned in the foreseeable future.
}
\references{
  \url{https://sourceforge.net/projects/liborigin/}: the \code{liborigin}
  implementation used by this package. Seems to be maintained as of 2019.
}
\keyword{package}
\examples{
  \dontrun{x <- read.opj('data.opj')}
}
