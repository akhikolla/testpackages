\docType{methods}
\name{setInputDefinitions,HybridSequenceClassifier-method}
\alias{setInputDefinitions,HybridSequenceClassifier-method}
\title{Sequence Detector Method: redefine all input definitions}
\usage{
\S4method{setInputDefinitions}{HybridSequenceClassifier}(fields, timestamp_start_field, timestamp_finish_field, 
context_field=NULL, preclassifier=NULL, preprocessor=NULL, pattern_field=NULL)
}
\arguments{
\item{fields}{(vector, character) - A vector of all relevant consolidated data stream fields.}
\item{timestamp_start_field}{(character) - A name of the field having starting time point values.}
\item{timestamp_finish_field}{(character) - A name of the field having finishing time point values.}
\item{context_field}{(character) - A name of the context identifier field (key field). If NULL, then \emph{.key} field is used for retrieving context identifier values.}
\item{preclassifier}{(\link{HSC_PC}) - A pre-classifier object. If NULL, the Sequence Detector creates new \emph{HSC_PC_None} pre-classifier, which means that the input consolidated data stream must have \emph{.clazz} field for retrieving classification values (input symbols in the underlying ETTs).}
\item{preprocessor}{(\link{HSC_PP}) - A pre-processing object. If NULL, the Sequence Detector creates new \emph{HSC_PP} pre-processor having the same fields as define in the \emph{fields} parameter, and ordering timestamp field as defined in \emph{timestamp_start_field}.}
\item{pattern_field}{(character) - A name of the field having output symbol values, i.e., relational ETT classification output.}
}
\description{
A method for redefining the Sequence Detector input parameters. This method is useful when we want to reuse an existing Sequence Detector for a different set of input data streams. Based on the ETT definition, after the pre-processing and pre-classification stages we need to have a consolidated data frame that comprises context identifier, sequence fields (timestamps or incremental value) and classification values (an input symbol). Not everything can be redefined and needs to be left as defined at the time of Sequence Detector instantiation, such as decay descriptors.
}

