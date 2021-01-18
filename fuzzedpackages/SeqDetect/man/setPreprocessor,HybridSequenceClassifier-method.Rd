\docType{methods}
\name{setPreprocessor,HybridSequenceClassifier-method}
\alias{setPreprocessor,HybridSequenceClassifier-method}
\title{Sequence Detector Method: re-set the pre-processor object}
\usage{
\S4method{setPreprocessor}{HybridSequenceClassifier}(preprocessor)
}
\arguments{
\item{preprocessor}{(\code{\link{HSC_PP}}) - New pre-processor object.}
}
\description{
Sequence Detector method for re-setting the pre-processor object. This might be desirable when we want to use already existing Sequence Detector for new input data streams, having different structure.
}
\seealso{
\code{\link{HybridSequenceClassifier}}
}
