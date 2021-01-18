\name{describe}
\Rdversion{1.1}
\alias{describe}
\title{Create descriptors to mutexes and attach}
\description{
The \code{describe} function returns information that is needed to 
``connect'' to a mutex from another process.  This connection is performed
by the \code{attach.mutex} function.}
\usage{
describe(x)
}
\arguments{
  \item{x}{a \code{boost.mutex} object}
}
\value{
The \code{describe} function returns a \code{boost.mutex.descriptor} object.
}
\author{Michael J. Kane <bigmemoryauthors@gmail.com>}
\examples{
  m = boost.mutex()
  mm = attach.mutex(describe(m))
  # Now, both m and mm specify the same mutex.
  rm(m)
  rm(mm)
  gc()
}
\keyword{programming}
\keyword{misc}
