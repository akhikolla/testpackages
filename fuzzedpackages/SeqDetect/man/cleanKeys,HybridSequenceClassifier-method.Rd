\docType{methods}
\name{cleanKeys,HybridSequenceClassifier-method}
\alias{cleanKeys,HybridSequenceClassifier-method}
\title{Sequence Detector Method: clean keys and tokens in machines (ETTs)}
\usage{
\S4method{cleanKeys}{HybridSequenceClassifier}(machine_id=NULL)
}
\arguments{
\item{machine_id}{(character) - An identifier of the machine (ETT) whose token needs to be removed. If NULL, all machines tokens are removed.}
}
\description{
Sequence Detector method for removing tokens.
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
hsc$cleanKeys()
}
