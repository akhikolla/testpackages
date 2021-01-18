\docType{methods}
\name{setPreclassifier,HybridSequenceClassifier-method}
\alias{setPreclassifier,HybridSequenceClassifier-method}
\title{Sequence Detector Method: re-set the pre-classifier object}
\usage{
\S4method{setPreclassifier}{HybridSequenceClassifier}(preclassifier)
}
\arguments{
\item{preclassifier}{(\code{\link{HSC_PC}}) - New pre-classifier object.}
}
\description{
Sequence Detector method for re-setting the pre-classifier object. This might be desirable when we want to use already existing Sequence Detector for new input data streams, having different structure.
}
\seealso{
\code{\link{HybridSequenceClassifier}}
}
