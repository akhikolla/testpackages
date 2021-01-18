\name{tennis}
\alias{tennis}
\alias{doubles}
\alias{doubles_noghost}
\alias{tennis_noghost}
\alias{doubles_ghost}
\alias{tennis_ghost}
\docType{data}
\title{Match outcomes from repeated doubles tennis matches}
\description{Match outcomes from repeated doubles tennis matches}
\usage{data(tennis)}
\format{
A hyper2 object corresponding to the match outcomes listed below.
}
\details{

There are four players, \eqn{p_1}{p1} to \eqn{p_4}{p4}.  These players
play doubles tennis matches with the following results:

\tabular{cc}{
  match \tab score\cr
  \eqn{\lbrace p_1,p_2\rbrace}{{p1,p2}} vs \eqn{\lbrace p_3,p_4\rbrace}{{p3,p4}} \tab 9-2\cr
  \eqn{\lbrace p_1,p_3\rbrace}{{p1,p3}} vs \eqn{\lbrace p_2,p_4\rbrace}{{p2,p4}} \tab 4-4\cr
  \eqn{\lbrace p_1,p_4\rbrace}{{p1,p4}} vs \eqn{\lbrace p_2,p_3\rbrace}{{p2,p3}} \tab 6-7\cr
  \eqn{\lbrace p_1\rbrace}{{p1}} vs \eqn{\lbrace p_3\rbrace}{{p3}} \tab 10-14\cr
  \eqn{\lbrace p_2\rbrace}{{p2}} vs \eqn{\lbrace p_3\rbrace}{{p3}} \tab 12-14\cr
  \eqn{\lbrace p_1\rbrace}{{p1}} vs \eqn{\lbrace p_4\rbrace}{{p4}} \tab 10-14\cr
  \eqn{\lbrace p_2\rbrace}{{p2}} vs \eqn{\lbrace p_4\rbrace}{{p4}} \tab 11-10\cr
  \eqn{\lbrace p_3\rbrace}{{p3}} vs \eqn{\lbrace p_4\rbrace}{{p4}} \tab 13-13\cr
}

It is suspected that \eqn{p_1}{p1} and \eqn{p_2}{p2} have some form of
team cohesion and play better when paired than when either solo or with
other players.  As the scores show, each player and, apart from p1-p2,
each doubles partnership, is of approximately the same strength.

Dataset \code{tennis} gives the appropriate likelihood function for the
players' strengths; and dataset \code{tennis_ghost} gives the
appropriate likelihood function if the extra strength due to team
cohesion of \eqn{\lbrace p_1,p_2\rbrace}{{p1,p2}} is represented by a
ghost player.
}
\source{Doubles tennis matches at \acronym{NOCS}, Jan-May 2008}
\references{
  Robin K. S. Hankin (2010).
  \dQuote{A Generalization of the Dirichlet Distribution},
  \emph{Journal of Statistical Software}, 33(11), 1-18,
  \url{http://www.jstatsoft.org/v33/i11/}
}
\examples{
data(tennis)
dotchart(maxp(tennis))
eigen(maxp(tennis,give=TRUE,hessian=TRUE)$hessian,TRUE,TRUE)$values

maxp(tennis_ghost)

}
\keyword{datasets}
