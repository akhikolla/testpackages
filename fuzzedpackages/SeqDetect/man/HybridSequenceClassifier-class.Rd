\docType{class}
\name{HybridSequenceClassifier-class}
\alias{HybridSequenceClassifier-class}
\alias{HybridSequenceClassifier}
\title{Sequence Detector}
\section{Fields}{
\describe{
\item{fields}{(vector, character) - A vector of all relevant consolidated data stream fields.}
\item{timestamp_start_field}{(character) - A name of the field having starting time point values.}
\item{timestamp_finish_field}{(character) - A name of the field having finishing time point values.}
\item{context_field}{(character) - A name of the context identifier field (key field). If NULL, then \emph{.key} field is used for retrieving context identifier values.}
\item{preclassifier}{(\link{HSC_PC}) - A pre-classifier object. If NULL, the Sequence Detector creates new \emph{HSC_PC_None} pre-classifier, which means that the input consolidated data stream must have \emph{.clazz} field for retrieving classification values (input symbols in the underlying ETTs).}
\item{preprocessor}{(\link{HSC_PP}) - A pre-processing object. If NULL, the Sequence Detector creates new \emph{HSC_PP} pre-processor having the same fields as define in the \emph{fields} parameter, and ordering timestamp field as defined in \emph{timestamp_start_field}.}
\item{decay_descriptors}{(list) - A list of decay descriptors. If NULL, token decay machanism is not used. Descriptor structure can be seen in vignettes.}
\item{pattern_field}{(character) - A name of the field having output symbol values, i.e., relational ETT classification output.}
\item{time_series_sequence_stats}{(logical) - If TRUE, ETTs are instructed to create sequence statistics. This is used whe having input time-series data streams. If FALSE, the sequence statistics are not created.}
\item{reuse_states}{(logical) - The parameter defined in [1]. ETTs are created so that each ETT have a state that represents each input symbol.}
\item{parallel_execution}{(logical) - Force parallel execution of ETTs in the Sequence Detector object. Useful when we expect higher number of ETTs in the same Sequence Detector.}
}}
\section{Methods}{
\describe{
\item{\code{cleanKeys(machine_id=NULL)}}{
  Sequence Detector method for removing tokens and keys\cr\code{\link{cleanKeys,HybridSequenceClassifier-method}}
}
\item{\code{clone()}}{
  Sequence Detector method for cloning\cr\code{\link{clone,HybridSequenceClassifier-method}}
}
\item{\code{compressMachines(ratio=0.5)}}{
  Sequence Detector method for compressing the underlying set of ETTs\cr\code{\link{compressMachines,HybridSequenceClassifier-method}}
}
\item{\code{getMachineIdentifiers()}}{
  Sequence Detector method for retrieving identifiers for the underlying set of ETTs\cr\code{\link{getMachineIdentifiers,HybridSequenceClassifier-method}}
}
\item{\code{induceSubmachine(threshold, isolate=FALSE)}}{
  Sequence Detector method for performing statistical projections on the underlying set of ETTs\cr\code{\link{induceSubmachine,HybridSequenceClassifier-method}}
}
\item{\code{mergeMachines()}}{
  Sequence Detector method for merging the underlying set of ETTs\cr\code{\link{mergeMachines,HybridSequenceClassifier-method}}
}
\item{\code{plotMachines(machine_id=NULL)}}{
  Sequence Detector method for plotting the underlying set of ETTs\cr\code{\link{plotMachines,HybridSequenceClassifier-method}}
}
\item{\code{printMachines(machine_id=NULL, state=NULL, print_cache=TRUE, print_keys=TRUE)}}{
  Sequence Detector method for printing the underlying set of ETTs to the R console\cr\code{\link{printMachines,HybridSequenceClassifier-method}}
}
\item{\code{process(streams, learn=TRUE, give_explain=TRUE, threshold=NULL, debug=FALSE, out_filename=NULL, ...)}}{
  Sequence Detector method for processing an input streams slice\cr\code{\link{process,HybridSequenceClassifier-method}}
}
\item{\code{serialize()}}{
  Sequence Detector method for serializing the underlying set of ETTs definitions\cr\code{\link{serialize,HybridSequenceClassifier-method}}
}
\item{\code{serializeToList()}}{
  Sequence Detector method for serializing the underlying set of ETTs definitions to the list\cr\code{\link{serializeToList,HybridSequenceClassifier-method}}
}
\item{\code{setOutputPattern(states=c(), transitions=c(), pattern, machine_id=NULL)}}{
  Sequence Detector method for setting the output alphabet to the underlying set of ETTs\cr\code{\link{setOutputPattern,HybridSequenceClassifier-method}}
}
\item{\code{setPreprocessor(preprocessor)}}{
  Sequence Detector method for setting the pre-processor\cr\code{\link{setPreprocessor,HybridSequenceClassifier-method}}
}
\item{\code{setPreclassifier(preclassifier)}}{
  Sequence Detector method for setting the pre-classifier\cr\code{\link{setPreclassifier,HybridSequenceClassifier-method}}
}
\item{\code{setInputDefinitions(fields, timestamp_start_field, timestamp_finish_field, context_field=NULL, preclassifier=NULL, preprocessor=NULL, pattern_field=NULL)}}{
  Sequence Detector method for redefining the input definitions\cr\code{\link{setInputDefinitions,HybridSequenceClassifier-method}}
}
}}
\description{
The Sequence Detector class.
}
\details{
Instantiates a Sequence Detector object. Constructor takes a number of parameters that define pre-processing
and pre-classification stages, as well as the structure of the input consolidated data stream.
These stages can be redefined again later using \code{\link{setInputDefinitions,HybridSequenceClassifier-method}} method.
See the SeqDetect vignette for examples.
}
\references{
[1] D. Krleža, B. Vrdoljak, and M. Brčić, Latent Process Discovery using Evolving Tokenized Transducer, \emph{IEEE Access}, vol. 7, pp. 169657 - 169676, Dec. 2019
}

