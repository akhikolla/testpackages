\docType{methods}
\name{compressMachines,HybridSequenceClassifier-method}
\alias{compressMachines,HybridSequenceClassifier-method}
\title{Sequence Detector Method: compress machines (ETTs)}
\usage{
\S4method{compressMachines}{HybridSequenceClassifier}(ratio=0.5)
}
\arguments{
\item{ratio}{(numeric) - A minimal isomorphic overlap between ETTs to be eligible for compression. Using this parameter too low (e.g. <0.5) might lead to overfragmentation of ETTs.}
}
\description{
Sequence Detector method for compressing machines by isolating common isomorphic sub-structured into child ETTs. See the SeqDetect vignette for details and examples.
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
hsc$process(input_streams)
hsc$compressMachines()
}
