\name{formula1}
\alias{formula1}
\alias{formula_1}
\alias{formula_one}
\alias{F1}
\alias{F1_2014}
\alias{F1_2015}
\alias{F1_2016}
\alias{F1_2017}
\docType{data}
\title{Formula 1 dataset}
\description{Race results from 2017 Formula One World Championship}
\usage{data(formula1)}
\format{A \code{hyper2} object that gives a likelihood function}
\details{

Object \code{formula1} is a \code{hyper2} object that gives a likelihood
function for the strengths of the competitors of the 2017 Formula One
World Championship.  It is created from \code{inst/formula1.txt} by the
code in \code{inst/formula1.R}, which is heavily documented.

}
\references{

  \dQuote{Wikipedia contributors}, \emph{2017 Formula One World
    Championship---{Wikipedia}{,} The Free Encyclopedia}, 2018.
  \url{https://en.wikipedia.org/w/index.php?title=2017_Formula_One_World_Championship&oldid=839923210}
  [Online; accessed 14-May-2018]
  }

\examples{
data(formula1)
dotchart(maxp(formula1))

}
