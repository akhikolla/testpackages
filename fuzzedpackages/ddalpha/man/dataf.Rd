\name{dataf}
\alias{dataf}
\docType{data}
\title{
Converts data from fdata class to the functional class.
}
\description{
\pkg{fda.usc} contains a handy function \code{\link[fda.usc:fdata]{fdata}} that converts varios types of functional data to the \code{fdata} class. 
To use these data in \code{\link{ddalphaf.train}} it must first be converted with \code{dataf}.

The function may be used either to convert a fdata object that contains multiple classes, or to convert multiple fdata objects, each of which contains one class.

Note that \code{fdata$fdata2d = TRUE} is not supported.
}
\usage{

dataf(fdatas, labels) 

}
\arguments{
  \item{fdatas}{ 
  an \code{fdata} object with curves belonging to multiple classes, or a list of \code{fdata} objects, each of which contains curves of the same class
}
  \item{labels}{
  a list of labels of the functional observations.
  If \code{fdatas} is a single \code{fdata} object, the list contains labels for each curve. 
  If \code{fdatas} is a list of \code{fdata} objects, the list labels for each of these \code{fdata} objects. 
}
}
\format{
  The functional data as a data structure (see \code{\link{dataf.*}}).
  \describe{
    \item{\code{dataf}}{
      The functional data as a list of objects. Each object is characterized by two coordinates 
      \describe{
        \item{\code{args}}{The arguments vector}
        \item{\code{vals}}{The values vector}
       }
    }
    \item{\code{labels}}{The classes of the objects}
  }
}
\seealso{
\code{\link{dataf.*}} for the functional data format.


 \code{\link{ddalphaf.train}} to train the functional DD\eqn{\alpha}-classifier
 
 \code{\link{compclassf.train}} to train the functional componentwise classifier


\code{\link{plot.functional}} for building plots of functional data
}
\examples{
\dontrun{
library(fda.usc)
data(phoneme)

# 1. convert a fdata object that contains multiple classes.
#    labels are defined for each curve
converted = dataf(phoneme$learn, phoneme$classlearn)
plot.functional(converted)

# 2. convert multiple fdata objects, each of which contains one class
#    the same label is applied to all curves of each fdata object
converted = dataf(list(phoneme$learn, phoneme$test), c("1 red", "2 blue"))
converted$name = "Phoneme learn (red) and test (blue)"
plot.functional(converted)
}
}
\keyword{fdata}
\keyword{functional}
