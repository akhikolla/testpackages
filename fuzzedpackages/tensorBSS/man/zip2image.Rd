\name{zip2image}
\alias{zip2image}
\title{ function to convert row of zip file to format used by image() }
\description{
  This .RD-file and the corresponding function are originally from the R-package ElemStatLearn which has now been removed from CRAN.

  This is a utility function converting zip.train/zip.test data 
  to format useful for plotting with the function \code{\link{image}}.
}
\usage{
zip2image(zip, line)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{zip}{ \code{\link{zip.train}} or \code{\link{zip.test}}.}
  \item{line}{ row of matrix to take }
}
\value{
  16 x 16 matrix suitable as argument for \code{\link{image}}.
}
\author{ Kjetil Halvorsen }
\references{
\cite{Kjetil B Halvorsen (package maintainer) (2019), R-package ElemStatLearn: Data Sets, Functions and Examples from the Book: "The Elements
of Statistical Learning, Data Mining, Inference, and Prediction" by Trevor Hastie, Robert Tibshirani and Jerome Friedman}
}
\examples{
## See example section of help file for zip.train
}
\keyword{ dplot }% 
