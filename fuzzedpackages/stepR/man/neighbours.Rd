\name{neighbours}
\alias{neighbours}
\alias{neighbors}

\title{Neighbouring integers}

\description{
Find integers within some radius of the given ones.
}

\usage{
neighbours(k, x = 1:max(k), r = 0)
}

\arguments{
\item{k}{integers within whose neighbourhood to look}
\item{x}{allowed integers}
\item{r}{radius within which to look}
}

\value{
Returns those integers in \code{x} which are at most \code{r} from some integer in \code{k}, i.e. the intersection of \code{x} with the union of the balls of radius \code{r} centred at the values of \code{k}. The return values are unique and sorted.
}

\seealso{\code{\link{is.element}}, \code{\link{match}}, \code{\link{findInterval}}, \code{\link{stepcand}}}

\examples{
neighbours(c(10, 0, 5), r = 1)
neighbours(c(10, 0, 5), 0:15, r = 1)
}

\keyword{nonparametric}
