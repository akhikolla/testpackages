\docType{methods}
\name{clone,HybridSequenceClassifier-method}
\alias{clone,HybridSequenceClassifier-method}
\title{Sequence Detector Method: clone the Sequence Detector object}
\usage{
\S4method{clone}{HybridSequenceClassifier}()
}
\description{
Sequence Detector method for cloning. Clones the Sequence Detector object and all its ETTs.
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

tt <- data.frame(product=c("P672","P113","P983","P23872","P5","P672","P2982","P983","P672",
                           "P991","P983","P113","P2982","P344"),
                 sales=c(2,11,12,98,8,18,298,16,24,25,18,16,43,101),alert=NA)
test_streams <- list(stream=tt)
hsc2 <- hsc$clone()
hsc2$process(test_streams,learn=FALSE)
}
