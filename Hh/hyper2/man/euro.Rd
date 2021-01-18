\name{euro}
\alias{euro}
\alias{euro2009}
\alias{eurovision}
\alias{Eurovision}
\alias{Eurovision2009}
\alias{Eurovision_song_contest}
\docType{data}
\title{Eurovision Song contest dataset}
\description{Voting patterns from Eurovision 2009}
\usage{data(euro)}
\format{
A \code{hyper2} object that gives a likelihood
function.
}
\details{

Object \code{euro2009} is a \code{hyper2} object that gives a likelihood
function for the skills of the 18 competitor countries in semi-final 1
of the 2009 Eurovision Song contest.  It is created by the
self-contained code in \code{inst/euro.R}, which is heavily documented.

The reason for choosing this particular dataset is that Pat Altham
(Statistical Laboratory, Cambridge) considered it with a view to
discover similarities between voters.  In the current analysis, the
likelihood function \code{euro2009} assumes their independence.

}
\references{
  \itemize{
    \item 
  Wikipedia contributors, \dQuote{Eurovision Song Contest
    2009---{Wikipedia}{,} The Free Encyclopedia}, 2018,
    \url{https://en.wikipedia.org/w/index.php?title=Eurovision_Song_Contest_2009&oldid=838723921}  [Online; accessed 13-May-2018].
    \item P. M. E. Altham, personal communication
  }
}
\examples{
data(euro)

dotchart(maxp(euro2009))

x1 <- loglik(euro2009,indep(maxp  (euro2009)))
x2 <- loglik(euro2009,indep(equalp(euro2009)))

## if the competitors have equal strength, the asymptotic distribution
##   of 2*(x1-x2) is chi-square with 17 degrees of freedom:

pchisq(2*(x1-x2),df=size(euro2009)-1,lower.tail=FALSE)
}
