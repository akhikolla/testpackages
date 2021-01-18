\name{rugby}
\Rdversion{1.1}
\alias{rugby}
\alias{super_rugby}
\alias{rugby2016}
\alias{super2016}
\alias{super2016g}
\docType{data}
\title{
Super Rugby 2016 data
}
\description{
Data from 2016 Super Rugby season.  Two \code{hyper2} objects,
\code{super2016} and \code{super2016g}.  The second includes a
\dQuote{ghost} team whose strength helps the home team.
}
\usage{data(super_rugby)}
\references{
  \itemize{
    \item
  Wikipedia contributors, \dQuote{Super Rugby,} \emph{Wikipedia, The
  Free Encyclopedia,} \url{
  https://en.wikipedia.org/w/index.php?title=Super_Rugby&oldid=758848475}
(accessed January 7, 2017).
\item 
Wikipedia contributors, \dQuote{2016 Super Rugby season,}
\emph{Wikipedia, The Free Encyclopedia,}
\url{https://en.wikipedia.org/w/index.php?title=2016_Super_Rugby_season&oldid=758243138} (accessed January 4, 2017).
}
  
}
\seealso{\code{\link{loglik}}}
\examples{
data(super_rugby)
loglik(super2016,indep(maxp(super2016)))
loglik(super2016g,indep(maxp(super2016g)))

# differ by more than 2, ghost is signficant!
}
\keyword{datasets}
