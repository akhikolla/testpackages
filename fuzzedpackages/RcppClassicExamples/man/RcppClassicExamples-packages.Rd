\name{RcppClassicExamples-package}
\alias{RcppClassicExamples-package}
\alias{RcppClassicExamples}
\docType{package}
\title{Examples for the deprecated Rcpp R/C++ Interface library API}
\description{
  This package shows some simple examples for the use of the deprecated
  \emph{classic} API from the first implementation of \pkg{Rcpp}. 

  Note that the \pkg{RcppClassic} package has been deprecated since
  2010, all new development should use the \pkg{Rcpp} package instead.
}
\details{
  The \pkg{Rcpp} package provides a number of C++ classes that ease
  access to C++ from R. This comprises both passing parameters to
  functions, as well as returning results back from C++ to R.

  Two APIs are supported. The first is an older API which was first
  introduced mostly in 2006 and extended in 2008. This interface is used
  by a few other packages and will be supported going forward in the
  \link[RcppClassic]{RcppClassic-package} package.
  
  A second and newer API that was started in 2009 offers more
  functionality, see the \link[Rcpp]{Rcpp-package} package.

  The \pkg{RcppExamples} package provides some simple examples for use
  of \pkg{Rcpp}.
}
\seealso{
  See the \link[RcppExamples]{RcppExamples-package} for examples of the
  recommended \pkg{Rcpp} API and \link[Rcpp]{Rcpp-package} for
  documentation on the recommended API to extend R with C++ code, while
  the deprecated \link[RcppClassic]{RcppClassic-package} documents the
  older, deprecated API.
}
\author{Dominick Samperi wrote the initial versions of Rcpp (and
  RcppTemplate) during 2005 and 2006.  Dirk Eddelbuettel made some
  additions, and became maintainer in 2008. Dirk Eddelbuettel and Romain
  Francois have been extending Rcpp since 2009.
}
\keyword{package}

