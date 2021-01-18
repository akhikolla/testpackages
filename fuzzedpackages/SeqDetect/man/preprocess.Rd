\name{preprocess}
\alias{preprocess}
\title{Pre-processing method}
\usage{
preprocess(x, streams, ...)
}
\arguments{
\item{x}{The pre-processor object.}
\item{streams}{A named list that comprises input data streams. Each input data stream is a data frame comprising fields declared while creating the \link{HSC_PP} object.}
\item{...}{An additional list of parameters that can be used by the pre-proccessor.}
}
\value{
Returns a list that comprises:
\itemize{
\item \emph{obj} - A returning pre-processor object. Passed in the subsequent invocation as \emph{x}.
\item \emph{res} - An output event stream. A resulting data frame representing the output event stream that is ordered according to the timestamp / sequence field and comprises all declared fields.
}
}
\description{
A method that all pre-processor classes need to implement. It is the code that aggregates and consolidates input data streams into one output event stream.
}
\details{
Input streams can be created as\preformatted{streams -> list(stream1=x1,stream2=x2,....)
}

where \emph{x1} is a data frame and \emph{stream1} is the name of the stream. All examnples can be seen in the SeqDetect vignette.
}
