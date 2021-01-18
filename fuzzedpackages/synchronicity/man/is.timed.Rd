\name{is.timed, timeout}
\alias{is.timed}
\alias{timeout}
\title{ Timeout operations for boost.mutex objects }
\description{
The \code{is.timed} function tells if a \code{boost.mutex} object has a timeout.
The \code{timeout} function tells how long a mutex will wait for a timeout.
}
\usage{
is.timed(m)
timeout(m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{ a \code{boost.mutex} object to get timeout information for }
}
\value{
\code{is.timed} returns \code{TRUE} if the object has a timeout and \code{FALSE}
otherwise.  If a timeout has been set \code{timeout} returns the number
of seconds a \code{boost.mutex} object will attempt to acquire a lock and
\code{NULL} otherwise.
}
\author{Michael J. Kane <bigmemoryauthors@gmail.com>}

\seealso{\code{\link{synchronicity}}}
\examples{
x = boost.mutex(timeout=5)
y = boost.mutex()
print(is.timed(x))
print(is.timed(y))
print(timeout(x))
print(timeout(y))
}
\keyword{programming}
\keyword{misc}
