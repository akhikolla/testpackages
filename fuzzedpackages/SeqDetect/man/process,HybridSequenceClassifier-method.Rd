\docType{methods}
\name{process,HybridSequenceClassifier-method}
\alias{process,HybridSequenceClassifier-method}
\title{Sequence Detector Method: processing input data streams}
\usage{
\S4method{process}{HybridSequenceClassifier}(streams,learn=TRUE,give_explain=TRUE,threshold=NULL,debug=FALSE,
out_filename=NULL, ...)
}
\arguments{
\item{streams}{(list, data.frame) - A named list that comprises input data streams. Each list element is a data frame that represents one input data stream.}
\item{learn}{(logical) - Are ETTs in the Sequence Detector extendable? If TRUE, the Sequence Detector learns new sequences from the supplied input data streams.}
\item{give_explain}{(logical) - Determines elements that will be returned by the method. If TRUE, output explanation and sequence statistical analysis will be returned as well.}
\item{threshold}{(integer) - Needed threshold for the pushing mechanism. Pushing will work only for transitions that are above the supplied threshold. If NULL, all transitions are taken in consideration.}
\item{debug}{(logical) - A switch for debug printout.}
\item{out_filename}{(character) - A filename where the consolidated data stream should be written. The written file is in the CSV format. If NULL, file writing is skipped.}
\item{...}{ - An additional list of parameters passed into pre-processor and pre-classifier.}
}
\description{
Sequence Detector method for processing of input data streams. See the SeqDetect vignette for proper usage and cases.
}
\value{
A list that comprises the following elements:
\itemize{
  \item{stream}{ - The consolidated stream.}
}
If \emph{give_explain} is TRUE then an additional element is:
\itemize{
  \item{explanation}{ - Actual and potential output symbols for each data item of the consolidated data stream.}
}
If \emph{give_explain} is TRUE and \emph{time_series_sequence_stats} is TRUE then an additional element is:
\itemize{
  \item{sequences}{ - The complete sequence statistics for the input time-series data.}
}
}
\seealso{
\code{\link{HybridSequenceClassifier}}
}
\examples{
st <- data.frame(product=c("P45","P134","P45","P134","P134","P45","P134"),
                 sales=c(2,12,18,16,18,24,8),
                 alert=c(NA,NA,NA,NA,NA,"Alert P45","Alert P134"))
input_streams <- list(stream=st)
pp <- HSC_PP(c("product","sales"),"sequence_id",auto_id=TRUE)
pc <- HSC_PC_Binning(0,100,40,"sales")
hsc <- HybridSequenceClassifier(c("product","sales","sequence_id"),
       "sequence_id","sequence_id","product",pc,pp)
res <- hsc$process(input_streams)
message(res)
}
