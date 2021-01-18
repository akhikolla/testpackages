\name{volley}
\alias{volley}
\alias{volleyball}
\alias{volleyball_matrix}
\title{Team sports}
\description{
Converts a matrix of games into a likelihood function.
}
\usage{
volley(M)
}

\arguments{
  \item{M}{Matrix}
}
\details{

The canonical example is the volleyball dataset.  Each row is a
volleyball game; each column is a player.  An entry of 0 means
\dQuote{on the losing side}, an entry of 1 means \dQuote{on the winning
side}, and an entry of \code{NA} means \code{did not play}.
}

\author{Robin K. S. Hankin}
\seealso{\code{\link{saffy}}}
\examples{

data(volleyball)
jj <- volley(volleyball_matrix)
stopifnot(jj == volleyball)
}