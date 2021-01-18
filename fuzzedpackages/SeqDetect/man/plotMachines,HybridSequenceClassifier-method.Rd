\docType{methods}
\name{plotMachines,HybridSequenceClassifier-method}
\alias{plotMachines,HybridSequenceClassifier-method}
\title{Sequence Detector Method: plot machines (ETTs)}
\usage{
\S4method{plotMachines}{HybridSequenceClassifier}(machine_id=NULL)
}
\arguments{
\item{machine_id}{(character) - A machine identifier that needs to be plotted. If NULL, all machines are plotted.}
}
\description{
Sequence Detector method for plotting of machines in the Sequence Detector object. Plotting is following the output symbols of the states and transitions. For machines that don't have a small output alphabet could not be plotted fully and correctly.
}
\seealso{
\code{\link{HybridSequenceClassifier}}
}
\examples{
ldf1 <- data.frame(product=c("P1","P1","P1","P1"),sequence_id=c(1,3,5,7),
                   sales=c(5,76,123,1),alert=c(NA,NA,NA,"Alert P1"))
ldf2 <- data.frame(product=c("P2","P2","P2","P2"),sequence_id=c(2,4,6,8),
                   sales=c(21,76,123,42),alert=c(NA,NA,NA,"Alert P2"))
input_streams <- list(stream1=ldf1,stream2=ldf2)
pp <- HSC_PP(c("product","sales","alert","sequence_id"),"sequence_id")
pc <- HSC_PC_Attribute("sales")
hsc <- HybridSequenceClassifier(c("sequence_id","product","sales","alert"),
                                         "sequence_id","sequence_id",context_field="product",
                                         preclassifier=pc,preprocessor=pp,reuse_states=TRUE,
                                         pattern_field="alert")
hsc$process(input_streams,learn=TRUE)
hsc$cleanKeys()
hsc$mergeMachines()
hsc$plotMachines()
}
