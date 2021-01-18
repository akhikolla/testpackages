\name{classify}
\alias{classify}
\title{Pre-classifying method}
\usage{
classify(x, stream, ...)
}
\arguments{
  \item{x}{(\code{\link{HSC_PC}}) - A pre-classifier object}
  \item{stream}{(data.frame) - An input event stream}
  \item{...}{An additional list of parameters needed for the used pre-classifier.}
}
\value{
(data.frame) - An output, a consolidated stream. Each row in the output data.frame must have \emph{.clazz} field, containing the row classification value.
}
\description{
An abstract method that needs to be implemented by classes that derive \code{\link{HSC_PC}}.
It performs classification on the input event stream.
See the SeqDetect vignette for details on how to implement a \code{\link{HSC_PC}} derived class.
}
