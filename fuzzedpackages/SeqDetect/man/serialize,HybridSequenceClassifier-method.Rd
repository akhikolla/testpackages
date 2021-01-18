\docType{methods}
\name{serialize,HybridSequenceClassifier-method}
\alias{serialize,HybridSequenceClassifier-method}
\title{Sequence Detector Method: serialize the Sequence Detector object}
\usage{
\S4method{serialize}{HybridSequenceClassifier}()
}
\description{
Sequence Detector method for serializing. User needs to serialize the Sequence Detector object before saving. If not performed, Sequence Detector C++ part of the object is not saved properly, and cannot be restored later.
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
hsc$serialize()
#saveRDS(hsc,"test.RDS")
# Previous line is commented due to the CRAN checking policies
}
