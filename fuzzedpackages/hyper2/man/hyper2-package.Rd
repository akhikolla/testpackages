\name{hyper2-package}
\alias{hyper2-package}
\docType{package}
\title{
\packageTitle{hyper2}
}
\description{
\packageDescription{hyper2}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{hyper2}
\packageIndices{hyper2}
  A generalization of the Dirichlet distribution, using a more
  computationally efficient method than the \pkg{hyperdirichlet}
  package.  The software is designed for the analysis of order
  statistics and team games.
}
\author{
\packageAuthor{hyper2}

Maintainer: \packageMaintainer{hyper2}
}
\references{
  \itemize{
    \item 
    R. K. S. Hankin (2010). \dQuote{A Generalization of the Dirichlet Distribution},
    \emph{Journal of Statistical Software}, 33(11), 1-18,
    \url{http://www.jstatsoft.org/v33/i11/}
    \item R. K. S. Hankin 2017. \dQuote{Partial Rank Data with the hyper2
      Package: Likelihood Functions for Generalized Bradley-Terry Models}.
    \emph{The R Journal} 9:2, pages 429-439.
  }
}
\keyword{ package }
\examples{
H  <- hyper2(list(1,2,3,1:2,2:3),1:5)

data(chess)
maxp(chess)

}
