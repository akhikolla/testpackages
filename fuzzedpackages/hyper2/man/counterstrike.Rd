\name{counterstrike}
\alias{counterstrike}
\alias{counterstrike_likelihood}
\docType{data}
\title{Counterstrike}
\description{
A kill-by-kill analysis of a counterstrike game.
}
\usage{data(counterstrike)}
\details{


E-sports are a form of competition using video games.  E-sports are
becoming increasingly popular, with high-profile tournaments attracting
over 400 million viewers, and prize pools exceeding US$20m.

Counter Strike: Global Offensive (CS-GO) is a multiplayer first-person
shooter game in which two teams of five compete in an immersive virtual
reality combat environment.

CS-GO is distinguished by the ability to download detailed gamefiles in
which every aspect of an entire match is recorded, it being possible to
replay the match at will.

Statistical analysis of such gamefiles is extremely difficult, primarily
due to complex gameplay features such as cooperative teamwork,
within-team communication, and real-time strategic fluidity.

It is the task of the statistician to make robust inferences from such
complex datasets, and here I analyze an influential match
between \dQuote{FaZe Clan} and \dQuote{Cloud9}, two of the most
successful E-sports syndicates of all time, when they competed at Boston
2018.  
  
Dataset \code{counterstrike} is a loglikelihood function for the
strengths of ten counterstrike players.  It can be created by executing
the code in \code{inst/counterstrike.R}, which also shows how to create
synthetic datasets for the purposes of comparison.

The probability model is identical to that of \code{NBA}: when a player
kills (scores), this is taken to be a success of the whole team rather
than the shooter.

The \code{counterstrike} object is created in
\code{inst/counterstrike.R} (but is called \dQuote{\code{H}} there, to
avoid terrible confusion).  There is quite a lot of extra documentation
in that file.

Counterstrike dataset kindly supplied by Zachary Hankin.


}
\references{
  \itemize{
    \item \url{https://www.youtube.com/watch?v=XKWzlG4jDnI}
    \item \url{https://en.wikipedia.org/wiki/FaZe_Clan}
    \item \url{https://en.wikipedia.org/wiki/Cloud9}
  }
}
\examples{

data(counterstrike)
dotchart(maxp(counterstrike))

}
\keyword{datasets}
