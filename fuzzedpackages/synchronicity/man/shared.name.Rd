\name{shared.name}
\alias{shared.name}
\title{The name of a mutex's shared resource}
\description{
This function returns the shared resource associated with a
\code{boost.mutex} object.
}
\usage{
shared.name(m)
}
\arguments{
  \item{m}{ a \code{boost.mutex} object }
}
\value{
A string specifying the shared resource associated with the given 
\code{boost.mutex} object.
}
\author{Michael J. Kane <bigmemoryauthors@gmail.com>}
\seealso{\code{\link{synchronicity}}}
\examples{
x = boost.mutex()
print(shared.name(x))
}
\keyword{programming}
\keyword{misc}
