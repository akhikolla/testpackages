\name{get.best}
\alias{get.best}

\title{
Find best fitting regmed model from regmed.grid object.
}
\description{
Find best fitting regmed model from regmed.grid object using minimum
BIC to select model.
}
\usage{
get.best(fit.grid)
}
\arguments{
  \item{fit.grid}{ a regmed.grid object
}
}
\value{
\item{ fit }{best fit regmed object based on minimum BIC}
\item{ grid }{row out of grid.data of regmed.grid object corresponding to best fit}
}
\references{
Schaid, DJ, Sinnwell JP. (2020) Penalized Models for Analysis of
Multiple Mediators. To appear in Genetic Epidemiology
}
\author{
Dan Schaid, Greg Jenkins, Jason Sinnwell
}
\seealso{
\code{\link{regmed.grid}}
}

\keyword{ models }% use one of  RShowDoc("KEYWORDS")
