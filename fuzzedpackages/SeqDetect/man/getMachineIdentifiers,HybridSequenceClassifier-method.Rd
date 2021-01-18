\docType{methods}
\name{getMachineIdentifiers,HybridSequenceClassifier-method}
\alias{getMachineIdentifiers,HybridSequenceClassifier-method}
\title{Sequence Detector Method: retrieve machine identifiers (ETTs)}
\usage{
\S4method{getMachineIdentifiers}{HybridSequenceClassifier}()
}
\description{
Sequence Detector method for retrieving list of machine identifiers.
}
\value{
(list) A list of strings, represeting machine identifiers.
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
message(hsc$getMachineIdentifiers())
}
