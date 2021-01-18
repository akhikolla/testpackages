\name{length.hyper2}
\alias{length.hyper2}
\alias{length}
\title{Length method for hyper2 objects}
\description{
Length method for hyper2 objects, being the number of different brackets
in the expression}
\usage{
\method{length}{hyper2}(x)
}
\arguments{
  \item{x}{hyper2 object}
}
\author{Robin K. S. Hankin}
\examples{

data("oneill")
length(icons)
seq_along(icons)
}
