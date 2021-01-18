\docType{methods}
\name{serializeToList,HybridSequenceClassifier-method}
\alias{serializeToList,HybridSequenceClassifier-method}
\title{Sequence Detector Method: serialize and externalize Sequence Detector object}
\usage{
\S4method{serializeToList}{HybridSequenceClassifier}()
}
\value{
Returns a list that comprises all Sequence Detector details. 
}
\description{
Sequence Detector method for serializing to a list. The list can be saved, loaded and deserialized into a Sequence Detector object again using \code{\link{deserializeFromList}} function.
}
\seealso{
\code{\link{HybridSequenceClassifier}},
\code{\link{deserializeFromList}}
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
hsc_list <- hsc$serializeToList()
#saveRDS(hsc_list,"test_list.RDS")
# Previous line is commented due to the CRAN checking policies
}
