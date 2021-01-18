\name{getPostExpectedgAndShrinkage}
\alias{getPostExpectedg}
\alias{getPostExpectedShrinkage}
\title{Extract posterior expected g and shrinkage factor (g/(1+g)) from a model}
\description{
  For a valid BayesMfp object, extract posterior expected g and
  shrinkage factor g/(1+g)
}
\usage{
getPostExpectedg(x, design=getDesignMatrix(x), nObs = nrow(design),
dim = ncol(design))
getPostExpectedShrinkage(x, design=getDesignMatrix(x),
nObs =nrow(design), dim = ncol(design))
}
\arguments{
  \item{x}{a valid \code{\link{BayesMfp}}-Object of length 1 (otherwise only first
    element is recognized)}
  \item{design}{the (centered) design matrix}
  \item{nObs}{number of observations}
  \item{dim}{number of design covariates}  
}
\value{
  The posterior expected g or shrinkage factor g/(1+g) from the model.
}
\author{Daniel Saban\'es Bov\'e}
\keyword{utilities}
\keyword{internal}
