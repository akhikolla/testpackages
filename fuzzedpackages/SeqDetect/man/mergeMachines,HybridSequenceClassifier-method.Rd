\docType{methods}
\name{mergeMachines,HybridSequenceClassifier-method}
\alias{mergeMachines,HybridSequenceClassifier-method}
\title{Sequence Detector Method: merge machines (ETTs)}
\usage{
\S4method{mergeMachines}{HybridSequenceClassifier}()
}
\description{
Sequence Detector method for merging machines. See the SeqDetect vignette for details and examples.
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
hsc$printMachines()
}
