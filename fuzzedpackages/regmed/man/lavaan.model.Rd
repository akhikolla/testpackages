\name{lavaan.model}
\alias{lavaan.model}

\title{
Create a lavaan model
}
\description{
Create a lavaan model for input to lavan::sem()
}
\usage{
lavaan.model(y.name, x.name, med.name, medcov)
}
\arguments{
  \item{y.name}{
character string of name of y variable
}
  \item{x.name}{
character string of name of x variable
}
  \item{med.name}{
vector of character strings of names for mediators
}
  \item{medcov}{
matrix of mediation covariances (square dimension to correspond with
med.name). This matrix is fixed in the sem model, and typically from a regmed object
}
}
\value{
a character string that descibes the mediation model
}
\author{
Dan Schaid, Greg Jenkins, Jason Sinnwell
}
\examples{
  ## NOT RUN:
  ## See vignette: regmed_demo
}
\keyword{models }

