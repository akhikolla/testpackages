\name{NBA}
\alias{NBA}
\alias{NBA_likelihood}
\alias{basketball}
\docType{data}
\title{Basketball dataset}
\description{
A point-by-point analysis of a basketball game
}
\usage{data(NBA)}
\details{

Dataset \code{NBA} is a dataframe contains a point-by-point analysis of
a basketball match.  Each row corresponds to a point scored.  The first
column is the time of the score, the second is the number of points
scored, the third shows which team had possession at the start of play,
and the fourth shows which team scored.

The other columns show the players.  Table entries show whether or not
that particular player was on the pitch when the point was scored.

Likelihood function \code{NBA_likelihood} is a \code{hyper2} object that
gives the log-likelihood function for this dataset.  It may be created
by executing \code{basketball.R}, in the \code{inst/} directory.

Note that the two \dQuote{ghost} players represent the effect of having
possession.
}
\references{
  \itemize{
    \item \url{http://www.espn.com/nba/playbyplay?gameId=400954514}
    \item \url{https://www.nbafullhd.com/24828-2/}
    \item \url{http://www.NBAbase.com/}
  }
}
\examples{

data(NBA)
dotchart(maxp(NBA_likelihood))

}
\keyword{datasets}
