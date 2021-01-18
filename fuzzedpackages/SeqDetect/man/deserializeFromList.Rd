\name{deserializeFromList}
\alias{deserializeFromList}
\title{Sequence Detector Method: deserialize Sequence Detector object from external list}
\usage{
deserializeFromList(l)
}
\arguments{
\item{l}{(list) - A list containing a Sequence Detector details.}
}
\value{
(HybridSequenceClassifier) - Returns a deserialized Sequence Detector object.
}
\description{
Sequence Detector method for deserializing from a list.
}
\seealso{
\code{\link{HybridSequenceClassifier-class}},
\code{\link{serializeToList,HybridSequenceClassifier-method}}
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

hsc_list <- hsc$serializeToList()
saveRDS(hsc_list,"test_list.RDS")


new_hsc_list <- readRDS("test_list.RDS")
file.remove("test_list.RDS")
hsc2 <- deserializeFromList(new_hsc_list)
}
