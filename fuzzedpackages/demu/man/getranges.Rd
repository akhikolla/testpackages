\name{getranges}
\alias{getranges}
\title{Get variable ranges from a design matrix.}
\usage{
getranges(design)
}
\arguments{
\item{design}{An \eqn{n \times p} matrix of input settings}

}
\value{
A \eqn{p \times 2} matrix with the lower and upper bounds (rounded to nearest integer value) of all \eqn{p} variables in the design matrix.
}
\description{
\code{getranges()} is a helper function to get the lower/upper bounds of variables in a design matrix, used for rescaling the inputs to the \eqn{[0,1]} hypercube.
}
\examples{
library(demu)

design=matrix(runif(10,1,5),ncol=2,nrow=5)
getranges(design)
}
