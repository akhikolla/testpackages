\name{freealg-package}
\alias{freealg-package}
\docType{package}
\title{
\packageTitle{freealg}
}
\description{
\packageDescription{freealg}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{freealg}
\packageIndices{freealg}

}
\author{
\packageAuthor{freealg}

Maintainer: \packageMaintainer{freealg}
}
\keyword{ package }
\examples{

a <- as.freealg("x+xyx")
b <- as.freealg("4x +XyX")  # upper-case interpreted as inverse

a+b
stopifnot(a+b==b+a)   # should be TRUE

a*b ==b*a # FALSE; noncommutative algebra

as.freealg("1+X+xy")^3

rfalg()
rfalg()^2


}
