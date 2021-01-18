\name{tidy}
\alias{tidy}
\alias{keep}
\alias{discard}
\title{Keep or discard players}
\description{Functionality to keep or discard subsets of the players in
  a hyper2 object}
\usage{
tidy(H)
keep(H, keep, tidy=TRUE)
discard(H, discard, tidy=TRUE)
}
\arguments{
  \item{H}{A \code{hyper2} object}
  \item{tidy}{Boolean, with default \code{TRUE} meaning to return a
    mathematically identical, but tidied, likelihood function}
  \item{keep,discard}{Players to keep or discard.  May be character or integer}
}
\details{

  Function \code{tidy(H)} returns a hyper2 object mathematically
  identical to \code{H} but with unused players (that is, players that
  do not appear in any bracket) removed.

  Functions \code{keep()} and \code{discard()} will either keep or
  discard players specified in the second argument.

}
\author{Robin K. S. Hankin}
\note{Function \code{tidy()} is very, very, inelegant}
\examples{

data("oneill")
maxp(icons)
discard(icons,c("OA","WAIS"))

data("skating")
maxp(skating)[1:4]      # numbers work, keep the first four skaters
maxp(keep(skating,1:4))

}