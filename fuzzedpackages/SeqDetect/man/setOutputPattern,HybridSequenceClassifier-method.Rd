\docType{methods}
\name{setOutputPattern,HybridSequenceClassifier-method}
\alias{setOutputPattern,HybridSequenceClassifier-method}
\title{Sequence Detector Method: assign output symbols}
\usage{
\S4method{setOutputPattern}{HybridSequenceClassifier}(states=c(),transitions=c(),pattern,machine_id=NULL)
}
\arguments{
\item{states}{(vector,character) - A character vector that comprises state identifiers. The supplied symbol (output alphabet, pattern parameter) is assigned to these states.}
\item{transitions}{(vector,character) - A character vector that comprises transition identifiers. The supplied symbol (output alphabet, pattern parameter) is assigned to these transitions.}
\item{pattern}{(character) - An output symbol, an element of the output alphabet, that needs to be assigned to supplied states and transitions.}
\item{machine_id}{(character) - If defined, the output symbol assignment applies only to the machine having this identifier. If NULL, the output symbol assignment applies to all machines (ETTs) in this Sequence Detector object.}
}
\description{
Sequence Detector method for assigning output symbols to states and transitions. See the SeqDetect vignette for proper usage and cases.
}
\seealso{
\code{\link{HybridSequenceClassifier}}
}
